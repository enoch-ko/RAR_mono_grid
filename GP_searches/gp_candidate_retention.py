"""Online exact retention for finite GP searches; no astrophysical computation."""
from collections import Counter
import math
import json
import numpy as np

DEFAULT_THRESHOLDS = {'2sigma_one_sided': .5 * math.erfc(2.0 / math.sqrt(2.0)),
                      '5sigma_one_sided': .5 * math.erfc(5.0 / math.sqrt(2.0))}


class CandidateRetention:
    """Keep a baseline, threshold minima, maximum P, and exact Pareto frontier.

    No epsilon dominance: arbitrarily close probabilities can bracket a cut.
    Original candidate IDs are stable; saved record indices are newly dense.
    This saves best-found candidates, never certifies a global optimum.
    """
    def __init__(self, thresholds=None, failure_examples=8):
        self.thresholds = dict(DEFAULT_THRESHOLDS if thresholds is None else thresholds)
        if any(not math.isfinite(p) or not 0 < p <= 1 for p in self.thresholds.values()):
            raise ValueError('retention probabilities must be finite in (0,1]')
        self.total = self.valid = 0
        self.baseline = self.maximum = None
        self.winners = {name: None for name in self.thresholds}
        self.frontier = []
        self.reasons = Counter()
        self.examples = []
        self.example_limit = failure_examples

    def add(self, record, xi):
        row = dict(record, original_candidate_index=self.total)
        self.total += 1
        is_valid = bool(row.get('valid', False))
        if is_valid and (not all(math.isfinite(row.get(k, float('nan')))
                                for k in ['J_ref', 'conditional_logp', 'marginal_logp'])
                         or row['J_ref'] < 0 or row['conditional_logp'] > 0
                         or row['marginal_logp'] > 0):
            raise ValueError('a valid candidate has invalid finite J/logP scores')
        xi = np.asarray(xi, dtype=float)
        if xi.ndim != 1 or not np.all(np.isfinite(xi)):
            raise ValueError('candidate profile must be a finite vector')
        candidate = (row, xi.copy())
        if row.get('source') == 'SPARC_constant_0.5':
            if self.baseline is not None:
                raise ValueError('duplicate baseline')
            self.baseline = candidate
        if not is_valid:
            reason = str(row.get('reason', 'unspecified'))
            tag = (reason if reason in ['signed_gbar_below_floor', 'nonfinite_signed_gbar']
                   else 'probability_warning' if reason.startswith('probability_warning')
                   else 'numerical_or_other_failure')
            self.reasons[tag] += 1
            if len(self.examples) < self.example_limit:
                self.examples.append(row)
            return
        self.valid += 1
        j, lp, index = row['J_ref'], row['marginal_logp'], row['original_candidate_index']
        for name, probability in self.thresholds.items():
            previous = self.winners[name]
            if lp >= math.log(probability) and (previous is None or
                    (j, index) < (previous[0]['J_ref'], previous[0]['original_candidate_index'])):
                self.winners[name] = candidate
        if self.maximum is None or (-lp, j, index) < (
                -self.maximum[0]['marginal_logp'], self.maximum[0]['J_ref'],
                self.maximum[0]['original_candidate_index']):
            self.maximum = candidate
        if any(r['J_ref'] <= j and r['marginal_logp'] >= lp for r, _ in self.frontier):
            return
        self.frontier = [(r, x) for r, x in self.frontier
                         if not (j <= r['J_ref'] and lp >= r['marginal_logp'])]
        self.frontier.append(candidate)
        self.frontier.sort(key=lambda item: (item[0]['J_ref'], -item[0]['marginal_logp']))

    def export(self):
        candidates = self.frontier + [self.baseline, self.maximum] + list(self.winners.values())
        unique = {r['original_candidate_index']: (r, x) for item in candidates
                  if item is not None for r, x in [item]}
        original_ids = sorted(unique)
        dense = {old: i for i, old in enumerate(original_ids)}
        records = [unique[old][0] for old in original_ids]
        profiles = [unique[old][1] for old in original_ids]
        frontier = [dense[r['original_candidate_index']] for r, _ in self.frontier]
        winners = {name: {'probability': p, 'log_probability': math.log(p),
                          'status': 'not_reached_in_finite_search' if self.winners[name] is None else 'reached',
                          'record_index': None if self.winners[name] is None else
                          dense[self.winners[name][0]['original_candidate_index']]}
                   for name, p in self.thresholds.items()}
        summary = {'evaluated': self.total, 'valid': self.valid, 'invalid': self.total-self.valid,
                   'retained': len(records), 'discarded': self.total-len(records),
                   'failure_counts': dict(self.reasons), 'failure_examples': self.examples,
                   'failure_examples_limit': self.example_limit,
                   'frontier_size': len(frontier), 'threshold_winners': winners,
                   'baseline_record_index': None if self.baseline is None else
                   dense[self.baseline[0]['original_candidate_index']],
                   'maximum_probability_record_index': None if self.maximum is None else
                   dense[self.maximum[0]['original_candidate_index']]}
        return records, profiles, frontier, summary


class RetentionStorageLimitError(RuntimeError):
    """Storage safety stop, never a physical infeasibility or complete search."""


def limited_json(value, byte_limit, indent=None):
    text = json.dumps(value, indent=indent, allow_nan=False) + '\n'
    size = len(text.encode('utf-8'))
    if size > byte_limit:
        raise RetentionStorageLimitError(
            f'retained JSON needs {size} bytes, exceeding {byte_limit}; '
            'no frontier truncation; choose a reviewed larger limit or search budget')
    return text
