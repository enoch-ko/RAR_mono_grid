"""Isolated finite-signed-gbar rank policy, 2026-10-16.

Only in-memory bindings in a NEW process change. Frozen production files stay
untouched. Positive disc M/L is required; net gbar may be negative or zero.
This is an ordering diagnostic, not positive-magnitude/circular MI compatibility.
No abs, clipping, node removal, nuisance-weight renormalisation or probability=0.
"""
import copy
import hashlib
import math
from itertools import permutations, product
from pathlib import Path
import numpy as np
import common_isotonic_reference as common
import gp_rar_search_retained_corrected as gp

ROOT = Path(__file__).resolve().parent
POLICY_ID = "finite_signed_gbar_rank_v1"
PARENT_PINS = {
    "common_isotonic_reference.py": "ef291d8777536d3e286f4381d1e0ca2774e390c6bb5c2cbe90d3e4d14a5f5863",
    "common_reference_driver.py": "d24087a17d08d11a7dd7ff9454ab7b8fe6b0f28e972e34b3e953af0f35f49ee8",
    "gp_rar_search_retained_corrected.py": "f3b84e65c38cb4d6cc9619b7c71773d4698521aeba7351a03e3338cee1076d70",
}
POLICY = {
    "id": POLICY_ID, "gbar_domain": "all finite signed values, including zero",
    "disc_domain": "finite strictly positive M/L within declared knot bounds",
    "nuisance_domain": "all original positive-weight GH nodes retained",
    "net_acceleration_floor_applied": False,
    "legacy_floor_parameter": "ignored; not a physical or numerical cut",
    "interpretation": "signed-rank diagnostic, not positive-magnitude RAR compatibility",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_parents():
    for name, expected in PARENT_PINS.items():
        if sha(ROOT / name) != expected:
            raise ValueError("signed policy parent changed: " + name)


def evidence():
    return {"policy": copy.deepcopy(POLICY), "adapter_sha256": sha(__file__),
            "parent_sha256": dict(PARENT_PINS)}


def log_axis_mask(gbar):
    """DISPLAY ONLY. Never feed this mask to scoring or marginalisation."""
    gbar = np.asarray(gbar, float)
    return np.isfinite(gbar) & (gbar > 0)


def physical_interval(mean, direction, lower, upper, disk_response=None,
                      nuisance_terms=None, floor=None):
    """Same closed-form amplitude interval, with DISC bounds only.

    Last three parameters preserve the parent call interface; they do not
    constrain signed net acceleration. Force nonfinites are checked at scoring.
    """
    mean, direction = np.asarray(mean, float), np.asarray(direction, float)
    if not (np.isfinite([lower, upper]).all() and 0 < lower < .5 < upper):
        raise ValueError("invalid strictly positive disc bounds")
    if mean.ndim != 1 or mean.shape != direction.shape:
        raise ValueError("matching disc and direction vectors required")
    q0 = np.concatenate((mean - lower, upper - mean))
    q1 = np.concatenate((direction, -direction))
    lo, hi = 0., np.inf
    for intercept, slope in zip(q0, q1):
        if not np.isfinite(intercept) or not np.isfinite(slope):
            return None
        if slope > 0:
            lo = max(lo, -intercept / slope)
        elif slope < 0:
            hi = min(hi, -intercept / slope)
        elif intercept < 0:
            return None
    if not np.isfinite(hi) or hi < lo or hi < 0:
        return None
    return max(0., lo), hi


class SignedProbabilityEvaluator(gp.ProbabilityEvaluator):
    def evaluate(self, xi):
        xi = np.asarray(xi, float)
        if not np.isfinite(xi).all() or np.any(xi <= 0):
            return {"valid": False, "reason": "nonpositive_or_nonfinite_disc_ml"}
        signed = self.D @ xi
        all_gbar = self.nuisance_terms + signed
        if not np.isfinite(all_gbar).all():
            return {"valid": False, "reason": "nonfinite_signed_gbar"}
        # Preserve signed W=R*g and every original nuisance weight.
        disk_term = signed / (1e6 / (self.r * 3.086e19))
        original = gp.mono.log_monotonic_probability_for_gbar
        gp.mono.log_monotonic_probability_for_gbar = self._ordered_logp
        try:
            conditional, marginal = gp.mono.evaluate_disk_term(
                disk_term, self.vgas, self.vbul, self.r, self.vobs, self.errv,
                self.args, self.gas_nodes, self.gas_weights,
                self.bulge_nodes, self.bulge_weights)
        finally:
            gp.mono.log_monotonic_probability_for_gbar = original
        if not conditional[1] or not marginal[1]:
            return {"valid": False, "reason": "conditional=" + str(conditional[2]) +
                    "; marginal=" + str(marginal[2])}
        messages = list(conditional[3]) + list(marginal[3])
        if messages:
            return {"valid": False, "reason": "probability_warning: " + " | ".join(messages)}
        if not np.isfinite([conditional[0], marginal[0]]).all() or max(conditional[0], marginal[0]) > 0:
            return {"valid": False, "reason": "invalid_log_probability"}
        return {"valid": True, "conditional_logp": float(conditional[0]),
                "marginal_logp": float(marginal[0]),
                "minimum_nuisance_gbar_m_s2": float(np.min(all_gbar)),
                "nonpositive_nuisance_entries": int(np.count_nonzero(all_gbar <= 0)),
                "positive_net_acceleration_at_all_nodes": bool(np.all(all_gbar > 0)),
                "gbar_domain_policy": POLICY_ID}


class SignedFrozenBaseline(common.FrozenBaseline):
    """Parent scoring, ties and weighted expectation; no gbar sign/floor gate."""
    def __init__(self, radius, error, nominal_gbar, node_gbar, weights,
                 settings=common.Settings(), integral=None):
        self.r = np.array(radius, dtype=float, copy=True)
        self.e = np.array(error, dtype=float, copy=True)
        self.nominal = np.array(nominal_gbar, dtype=float, copy=True)
        self.nodes = np.array(node_gbar, dtype=float, copy=True)
        self.weights = np.array(weights, dtype=float, copy=True)
        self.settings = settings
        if (self.r.ndim != 1 or len(self.r) == 0 or self.e.shape != self.r.shape or
                self.nominal.shape != self.r.shape or self.nodes.ndim != 2 or
                self.nodes.shape[1:] != self.r.shape or self.weights.shape != (len(self.nodes),)):
            raise ValueError("invalid baseline shapes")
        if not all(np.isfinite(a).all() for a in (self.r, self.e, self.nominal, self.nodes, self.weights)):
            raise ValueError("nonfinite baseline")
        if np.any(self.r <= 0) or np.any(self.e <= 0):
            raise ValueError("invalid radius/error baseline")
        # Retain legacy settings validation, but floor_relative is NOT applied.
        if (settings.ngrid < 2 or settings.max_refinements < 1 or settings.max_tie_permutations < 1 or
                not np.isfinite([settings.zpadding, settings.log_tolerance, settings.floor_relative]).all() or
                min(settings.zpadding, settings.log_tolerance, settings.floor_relative) <= 0):
            raise ValueError("invalid numerical settings")
        if np.any(self.weights <= 0) or not np.isclose(np.sum(self.weights), 1., rtol=1e-12, atol=1e-15):
            raise ValueError("positive normalized nuisance weights required")
        # Deliberately do not drop nodes or renormalise supplied weights.
        self.orders = []
        for node in self.nodes:
            order = np.argsort(-node, kind="stable")
            groups = np.split(order, np.flatnonzero(np.diff(node[order]) != 0) + 1)
            count = math.prod(math.factorial(len(g)) for g in groups)
            if count > settings.max_tie_permutations:
                raise ValueError("tie cap exceeded; no stable-order fallback")
            self.orders.append(tuple(tuple(j for part in parts for j in part)
                for parts in product(*(permutations(g) for g in groups))))
        if integral is None:
            from probability_recursion import boundary_log_probability
            integral = boundary_log_probability
        self.integral = integral
        for a in (self.r, self.e, self.nominal, self.nodes, self.weights):
            a.setflags(write=False)


def signed_context(context):
    context = copy.deepcopy(context)
    context["signed_rank_policy"] = copy.deepcopy(POLICY)
    context.setdefault("invariants", {})["signed_rank_policy"] = copy.deepcopy(POLICY)
    return context


def install_gp():
    """Call only inside a dedicated new-process launcher, never a live worker."""
    verify_parents()
    if getattr(gp, "_signed_rank_policy_installed", False):
        raise RuntimeError("policy already installed")
    gp.physical_interval = physical_interval
    gp.ProbabilityEvaluator = SignedProbabilityEvaluator
    gp.SOURCE_FILES = (*gp.SOURCE_FILES, Path(__file__).name, "run_signed_rank_20261016.py")
    parent_run = gp.galaxy_run
    def run(galaxy, args, reference_path=None):
        if reference_path is not None:
            raise ValueError("legacy QP profile import not approved for signed-rank revision")
        args.signed_rank_policy = POLICY_ID
        args.net_acceleration_floor_applied = False
        result = parent_run(galaxy, args, None)
        verify_parents()
        result["signed_rank_policy"] = evidence()
        result["settings"]["signed_gbar_floor_relative"] = None
        result["physical_domain"]["policy"] = copy.deepcopy(POLICY)
        result["physical_domain"]["floor_m_s2"] = None
        result["diagnostic_only"] = True
        return result
    gp.galaxy_run = run
    gp._signed_rank_policy_installed = True


def install_reference():
    verify_parents()
    if getattr(common, "_signed_rank_policy_installed", False):
        raise RuntimeError("policy already installed")
    common.FrozenBaseline = SignedFrozenBaseline
    common.REFERENCE_ID = "common_baseline_isotonic_signed_rank_v1"
    parent_run, parent_save = common.run_reference, common.save_reference
    def run(velocity, scorer, *, B, seed, context, save_vectors=False):
        result = parent_run(velocity, scorer, B=B, seed=seed,
                            context=signed_context(context), save_vectors=save_vectors)
        result["signed_rank_policy"] = evidence()
        result["nonpositive_nuisance_entries"] = int(np.count_nonzero(scorer.nodes <= 0))
        result["nominal_log_axis_mask"] = log_axis_mask(scorer.nominal).tolist()
        return result
    def save(result, output, max_bytes=8*1024*1024):
        verify_parents()
        result.setdefault("driver_source_sha256", {})[Path(__file__).name] = sha(__file__)
        result["driver_source_sha256"]["run_signed_rank_20261016.py"] = sha(ROOT / "run_signed_rank_20261016.py")
        result.setdefault("build_plan", {})["signed_rank_policy"] = copy.deepcopy(POLICY)
        return parent_save(result, output, max_bytes=max_bytes)
    common.run_reference, common.save_reference = run, save
    common._signed_rank_policy_installed = True
