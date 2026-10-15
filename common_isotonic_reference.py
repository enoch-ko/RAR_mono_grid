"""Frozen SPARC-baseline empirical benchmark; not a calibrated significance test.
Standalone core only: no production launcher, candidate refits, or force solves.
Callers supply the constant-disc=.5 nominal and nuisance-node accelerations.
"""
from dataclasses import asdict, dataclass
from itertools import permutations, product
from pathlib import Path
import hashlib
import json
import math
import warnings
import numpy as np

PARSER_ID = "sparc_rotmod_header_v1_gas3_disk4"
REFERENCE_ID = "common_baseline_isotonic_v1"


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                    separators=(",", ":")).encode()).hexdigest()


def read_rotmod(path, radius_bounds):
    """Validate SPARC column names; never inherit the inverse loader's gas index."""
    path = Path(path)
    content = path.read_bytes()
    header = next((line.lstrip("# ").split() for line in content.decode().splitlines()
                   if line.startswith("#") and "Vobs" in line), None)
    expected = ["Rad", "Vobs", "errV", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]
    if header != expected:
        raise ValueError("unrecognized rotmod column header")
    raw = np.loadtxt(path, ndmin=2)
    lo, hi = map(float, radius_bounds)
    if raw.shape[1] != 8 or not (np.isfinite([lo, hi]).all() and 0 <= lo < hi):
        raise ValueError("invalid rotmod shape or photometric support")
    mask = np.isfinite(raw[:, 0]) & (raw[:, 0] > 0) & (raw[:, 0] >= lo) & (raw[:, 0] <= hi)
    retained = raw[mask]
    if len(retained) < 1 or not np.isfinite(retained[:, :6]).all() or np.any(retained[:, 2] <= 0):
        raise ValueError("invalid retained kinematic sample")
    return {name: retained[:, j].copy() for j, name in enumerate(expected)}, {
        "parser_id": PARSER_ID, "rotmod_sha256": hashlib.sha256(content).hexdigest(),
        "mask": mask.tolist(), "radius_bounds_kpc": [lo, hi]}


def isotonic_reference(gbar, radius, velocity, error):
    """Weighted nonnegative weak-order PAVA, pooling exact abscissa ties first."""
    x, r, v, e = [np.asarray(a, dtype=float) for a in (gbar, radius, velocity, error)]
    if x.ndim != 1 or len(x) == 0 or any(a.shape != x.shape for a in (r, v, e)):
        raise ValueError("nonempty matching vectors required")
    if not all(np.isfinite(a).all() for a in (x, r, v, e)) or np.any(r <= 0) or np.any(e <= 0):
        raise ValueError("finite inputs and positive radii/errors required")
    y, sigma = v / np.sqrt(r), e / np.sqrt(r)
    weights = (np.min(sigma) / sigma) ** 2
    if np.any(weights == 0):
        raise ValueError("weight dynamic range underflow")
    order = np.argsort(x, kind="stable")
    groups = np.split(order, np.flatnonzero(np.diff(x[order]) != 0) + 1)
    blocks = []
    for group in groups:
        blocks.append([list(group), float(np.sum(weights[group])),
                       float(np.dot(weights[group], y[group]))])
        while len(blocks) > 1 and blocks[-2][2] / blocks[-2][1] > blocks[-1][2] / blocks[-1][1]:
            b = blocks.pop(); a = blocks.pop()
            blocks.append([a[0] + b[0], a[1] + b[1], a[2] + b[2]])
    theta = np.empty_like(y)
    for indices, weight, weighted_y in blocks:
        theta[indices] = max(0.0, weighted_y / weight)
    return {"theta": theta.tolist(), "velocity_ref": (np.sqrt(r) * theta).tolist(),
            "nominal_gbar": x.tolist(), "ascending_order": order.tolist(),
            "exact_tie_groups": [g.tolist() for g in groups if len(g) > 1]}


@dataclass(frozen=True)
class Settings:
    ngrid: int = 400
    zpadding: float = 8.0
    log_tolerance: float = 2e-5
    max_refinements: int = 8
    max_tie_permutations: int = 720
    floor_relative: float = 1e-8


class FrozenBaseline:
    """Precompute node orderings once. Integral values live ONLY inside each call.

    node_gbar must already include the declared GH gas/bulge treatment and .5 disc.
    The floor is reference availability QC, not a new grid-candidate restriction.
    """
    def __init__(self, radius, error, nominal_gbar, node_gbar, weights,
                 settings=Settings(), integral=None):
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
        if np.any(self.r <= 0) or np.any(self.e <= 0) or np.any(self.nominal <= 0):
            raise ValueError("invalid radius/error/nominal baseline")
        if (settings.ngrid < 2 or settings.max_refinements < 1 or settings.max_tie_permutations < 1 or
                not np.isfinite([settings.zpadding, settings.log_tolerance, settings.floor_relative]).all() or
                min(settings.zpadding, settings.log_tolerance, settings.floor_relative) <= 0):
            raise ValueError("invalid numerical settings")
        if np.any(self.nodes <= 0) or np.any(self.nodes < settings.floor_relative * np.max(self.nominal)):
            raise ValueError("reference unavailable: nuisance-node acceleration floor")
        if np.any(self.weights <= 0) or not np.isclose(np.sum(self.weights), 1.0, rtol=1e-12, atol=1e-15):
            raise ValueError("positive normalized nuisance weights required")
        self.weights /= self.weights.sum()
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

    def score(self, velocity):
        from scipy.special import logsumexp
        v = np.asarray(velocity, dtype=float)
        if v.shape != self.r.shape or not np.isfinite(v).all():
            raise ValueError("invalid complete velocity vector")
        memo, values = {}, []
        for orders in self.orders:
            terms = []
            for order in orders:
                if order not in memo:
                    o = np.asarray(order)
                    with warnings.catch_warnings():
                        warnings.simplefilter("error")
                        logp = float(self.integral(v[o], self.e[o], self.r[o],
                            ngrid=self.settings.ngrid, zpadding=self.settings.zpadding,
                            log_tolerance=self.settings.log_tolerance,
                            max_refinements=self.settings.max_refinements))
                    if not np.isfinite(logp) or logp > 0:
                        raise ValueError("invalid conditional log probability")
                    memo[order] = logp
                terms.append(memo[order])
            values.append(float(logsumexp(terms) - np.log(len(terms))))
        from utils import log_weighted_average
        return float(log_weighted_average(np.asarray(values), self.weights))


class EmpiricalCDF:
    def __init__(self, log_scores):
        self.sorted = np.sort(np.array(log_scores, dtype=float, copy=True))
        if self.sorted.ndim != 1 or len(self.sorted) == 0 or not np.isfinite(self.sorted).all() or np.any(self.sorted > 0):
            raise ValueError("only complete finite backend log scores accepted")
        self.sorted.setflags(write=False)

    def evaluate(self, logp):
        if not np.isfinite(logp) or logp > 0:
            raise ValueError("finite candidate log probability required")
        n = len(self.sorted); k = int(np.searchsorted(self.sorted, logp, side="right"))
        q, z = k / n, 1.959963984540054
        den = 1 + z*z/n
        centre = (q + z*z/(2*n))/den
        half = z*np.sqrt(q*(1-q)/n + z*z/(4*n*n))/den
        return {"q_common_ref": q, "count_le": k, "B": n,
                "mc_resolution": 1/n, "conditional_mc_wilson95":
                [max(0.0, centre-half), min(1.0, centre+half)]}

    def cutoff(self, alpha):
        if not np.isfinite(alpha) or not 0 <= alpha <= 1:
            raise ValueError("alpha must be in [0,1]")
        if alpha == 0:
            return {"accept_all_valid_scores": True, "logp_cutoff": None}
        n = len(self.sorted); k = max(1, min(n, math.ceil(alpha*n)))
        while k < n and k/n < alpha:
            k += 1
        while k > 1 and (k-1)/n >= alpha:
            k -= 1
        return {"accept_all_valid_scores": False, "logp_cutoff": float(self.sorted[k-1])}


def run_reference(velocity, scorer, *, B, seed, context, save_vectors=False):
    """Any failed observation/mock disables the entire CDF; never reduce B."""
    required = {"galaxy", "parser_id", "input_sha256", "normalization_id", "nuisance_model"}
    if not required.issubset(context) or context["parser_id"] != PARSER_ID:
        raise ValueError("explicit compatible input/parser/normalization context required")
    if isinstance(B, bool) or not isinstance(B, (int, np.integer)) or B < 1:
        raise ValueError("B must be an explicit positive integer")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be an explicit nonnegative integer")
    fit = isotonic_reference(scorer.nominal, scorer.r, velocity, scorer.e)
    from utils import vel2acc
    fit["reference_acceleration_m_s2"] = vel2acc(np.asarray(fit["velocity_ref"])**2, scorer.r).tolist()
    rng = np.random.Generator(np.random.PCG64(int(seed)))
    initial = rng.bit_generator.state
    scores, failures, vectors = [], [], []
    try:
        observed = scorer.score(velocity)
    except TimeoutError:
        raise
    except Exception as exc:
        observed = None; failures.append({"stage": "observation", "reason": str(exc)})
    negative = 0
    for b in range(int(B)):
        draw = np.asarray(fit["velocity_ref"]) + scorer.e*rng.standard_normal(len(scorer.r))
        negative += int(np.any(draw < 0))
        if save_vectors:
            vectors.append(draw.tolist())
        try:
            scores.append(scorer.score(draw))
        except TimeoutError:
            raise
        except Exception as exc:
            scores.append(None); failures.append({"stage": "mock", "index": b, "reason": str(exc)})
    frozen_context = json.loads(json.dumps(context, allow_nan=False))
    frozen_context["probability_settings"] = asdict(scorer.settings)
    frozen_context["reference_id"] = REFERENCE_ID
    frozen_context["baseline_arrays_sha256"] = digest({"r": scorer.r.tolist(),
        "error": scorer.e.tolist(), "observed": np.asarray(velocity).tolist(),
        "nominal_gbar": scorer.nominal.tolist(), "node_gbar": scorer.nodes.tolist(),
        "weights": scorer.weights.tolist()})
    root = Path(__file__).resolve().parent
    frozen_context["source_sha256"] = {name: hashlib.sha256((root/name).read_bytes()).hexdigest()
        for name in ("common_isotonic_reference.py", "probability_recursion.py", "utils.py")}
    out = {"status": "unusable" if failures else "complete", "reference": fit,
           "context": frozen_context, "context_sha256": digest(frozen_context),
           "radii_kpc": scorer.r.tolist(), "quoted_error": scorer.e.tolist(),
           "observed_velocity": np.asarray(velocity).tolist(), "observed_logp": observed,
           "node_gbar": scorer.nodes.tolist(), "node_weights": scorer.weights.tolist(),
           "B": int(B), "seed": int(seed), "rng": "PCG64", "initial_rng_state": initial,
           "final_rng_state": rng.bit_generator.state, "mock_logp": scores,
           "negative_velocity_draws": negative, "failures": failures}
    if save_vectors:
        out["mock_velocity"] = vectors
    if not failures:
        out["sorted_mock_logp"] = sorted(scores)
        out["degenerate_point_mass"] = len(set(scores)) == 1
        out["display_underflow_count"] = sum(float(np.exp(x)) == 0 for x in scores)
    return out


def annotate(reference, candidate_logp, candidate_context):
    """Lookup only. Parser identity must match; input hashes alone are insufficient."""
    if reference["status"] != "complete":
        raise ValueError("failed/incomplete mock set cannot annotate candidates")
    if digest(reference["context"]) != reference["context_sha256"]:
        raise ValueError("reference context hash mismatch")
    if candidate_context != reference["context"]:
        raise ValueError("incompatible galaxy/input/parser/normalization/numerics")
    if candidate_context.get("parser_id") != PARSER_ID:
        raise ValueError("unverified input parser")
    if len(reference["mock_logp"]) != reference["B"] or reference.get("failures"):
        raise ValueError("incomplete mock denominator")
    return {"original_logp": float(candidate_logp),
            **EmpiricalCDF(reference["mock_logp"]).evaluate(candidate_logp)}


def prepare_sparc_baseline(data, disk_w_unit, settings=Settings()):
    """One precomputed Casertano W at M/L=1; no force solve during mock scoring.

    The caller must supply the declared normalized full-density disc solve,
    interpolated onto these retained radii; never substitute rotmod Vdisk.
    """
    from utils import lognormal_quadrature, vel2acc
    r = np.asarray(data["Rad"], dtype=float)
    w = np.asarray(disk_w_unit, dtype=float)
    if w.shape != r.shape or not np.isfinite(w).all():
        raise ValueError("matching finite unit-disc signed W required")
    gas = np.asarray(data["Vgas"]); bulge = np.asarray(data["Vbul"])
    gn, gw = lognormal_quadrature(1.0, .04, 9)
    bn, bw = lognormal_quadrature(.7, .1, 9)
    nominal = vel2acc(.5*w + gas*np.abs(gas) + .7*bulge**2, r)
    nodes = [vel2acc(.5*w + g*gas*np.abs(gas) + b*bulge**2, r)
             for g in gn for b in bn]
    return FrozenBaseline(r, data["errV"], nominal, nodes,
                          [g*b for g in gw for b in bw], settings)


def save_reference(result, output, max_bytes=8*1024*1024):
    """Exclusive new artifact, bounded retention; never overwrite old results."""
    payload = json.dumps(result, sort_keys=True, allow_nan=False).encode()
    if len(payload) > max_bytes:
        raise ValueError("reference artifact exceeds storage cap")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as stream:
        stream.write(payload)
    return hashlib.sha256(payload).hexdigest()
