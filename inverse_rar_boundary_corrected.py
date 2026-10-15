#!/usr/bin/env python
"""Isolated inverse diagnostic using boundary quadrature, leaving old jobs untouched.

python3 inverse_rar_boundary.py --galaxy UGC02916 --dry-run
Use the scheduler for full 9x9 nuisance evaluation; do not overwrite old results.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np
import probability_recursion as probability


def json_safe(value):
    """Preserve invalid flags/reasons while representing unavailable numbers as null."""
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(value):
        return None
    return value


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--galaxy', required=True)
    parser.add_argument('--boundary-log-tolerance', type=float, default=2e-5)
    parser.add_argument('--boundary-max-refinements', type=int, default=8)
    parser.add_argument('--dry-run', action='store_true')
    args, remaining = parser.parse_known_args()
    if args.boundary_log_tolerance <= 0 or not np.isfinite(args.boundary_log_tolerance):
        parser.error('positive finite boundary tolerance required')
    if args.boundary_max_refinements < 1:
        parser.error('at least one boundary refinement required')
    import inverse_rar_diagnostic_corrected as inverse
    original_analyze = inverse.analyze
    original_solver = probability.boundary_log_probability
    calls = []
    forwarded = ['--galaxy', args.galaxy] + remaining
    if not any(a == '--output' or a.startswith('--output=') for a in remaining):
        target = (Path(inverse.DEFAULT_PROJECT_ROOT)/'outputs'/'rotmod_columns_corrected_20261015'
                  /'qp'/f'{args.galaxy}_inverse.json')
        forwarded += ['--output', str(target)]
    settings = {'method': probability.METHOD,
                'log_tolerance': args.boundary_log_tolerance,
                'max_refinements': args.boundary_max_refinements,
                'resolution_meaning': 'starting neighbourhood/background points; doubled adaptively',
                'tail_convergence': 'must be checked independently of resolution refinement',
                'physics_cache_gate': 'unchanged disk-basis gate; old probability QC is not new-solver QC'}
    if args.dry_run:
        print(json.dumps({'forwarded_arguments': forwarded, 'probability_settings': settings}, indent=2))
        return

    def recorded_solver(*values, **options):
        options.update(log_tolerance=args.boundary_log_tolerance,
                       max_refinements=args.boundary_max_refinements, return_diagnostics=True)
        try:
            info = original_solver(*values, **options)
        except Exception as exc:
            calls.append({'valid': False, 'reason': f'{type(exc).__name__}: {exc}'})
            raise
        calls.append({'valid': True, 'log_probability': info['log_probability'],
                      'levels': info['levels'], 'zpadding': info['zpadding']})
        return info['log_probability']

    def analyzed(inverse_args):
        result = original_analyze(inverse_args)
        result['probability_numerics'] = dict(settings, ordered_solver_calls=calls)
        result['provenance']['boundary_entry_sha256'] = digest(__file__)
        result['provenance']['probability_recursion_sha256'] = digest(probability.__file__)
        result['provenance']['utils_sha256'] = digest(Path(inverse.DEFAULT_PROJECT_ROOT)/'utils.py')
        return json_safe(result)

    # Injection is confined to this NEW process. Existing batch processes and
    # on-disk inverse/mono/utils sources keep their original solver and hashes.
    probability.boundary_log_probability = recorded_solver
    inverse.mono.log_monotonic_probability_for_gbar = probability.boundary_log_probability_for_gbar
    inverse.analyze = analyzed
    sys.argv = [sys.argv[0]] + forwarded
    inverse.main()


if __name__ == '__main__':
    main()
