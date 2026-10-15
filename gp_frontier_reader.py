"""Bounded-memory reader of pretty-printed legacy GP frontiers and their records.

It verifies metadata and retained record syntax/counts, not every discarded
record's JSON syntax. It never alters files or computes forces/probabilities.
"""
import hashlib
import json
from pathlib import Path
import re


class UnsupportedFrontierLayout(ValueError):
    pass


def read_frontier_result(path, chunk_size=1024*1024):
    path = Path(path)
    initial = path.stat()
    signature = (initial.st_size, initial.st_mtime_ns)
    record_marker = b'\n  "records": ['
    frontier_marker = b'\n  "frontier": ['
    with path.open('rb') as stream:
        head = stream.read(1024*1024)
        if record_marker not in head:
            raise UnsupportedFrontierLayout('legacy root records marker unavailable')
        start = head.index(record_marker)
        record_start = start + len(record_marker)
        tail_size = min(initial.st_size, 1024*1024)
        while True:
            tail_offset = initial.st_size-tail_size
            stream.seek(tail_offset)
            tail = stream.read()
            frontier_start = tail.find(frontier_marker)
            if frontier_start >= 0:
                break
            if tail_size >= min(initial.st_size, 16*1024*1024):
                raise UnsupportedFrontierLayout('frontier metadata exceeds tail-reader budget')
            tail_size = min(initial.st_size, tail_size*2)
        closing = tail.rfind(b'\n  ],', 0, frontier_start)
        if closing < 0:
            raise UnsupportedFrontierLayout('records closing delimiter unavailable')
        record_end = tail_offset + closing
        if record_end < record_start:
            raise UnsupportedFrontierLayout('invalid record span')
        data = json.loads(head[:start] + b'\n  "records": [],' + tail[frontier_start:])
        if data.get('checkpoint_storage', {}).get('mode') == 'retained-only':
            raise UnsupportedFrontierLayout('retained mode should use the normal small JSON reader')
        original_frontier = [int(item['record_index']) for item in data['frontier']]
        wanted = set(original_frontier) | {0}
        found, index, buffer, offset = {}, 0, b'', 0
        digest = hashlib.sha256()
        stream.seek(0)
        pattern = re.compile(rb'\n    \{')
        while True:
            block = stream.read(chunk_size)
            if not block:
                break
            digest.update(block)
            lo, hi = max(record_start-offset, 0), min(record_end-offset, len(block))
            if hi > lo:
                buffer += block[lo:hi]
                positions = [match.start() for match in pattern.finditer(buffer)]
                for left, right in zip(positions, positions[1:]):
                    if index in wanted:
                        found[index] = json.loads(buffer[left:right].strip().rstrip(b',').strip())
                    index += 1
                if positions:
                    buffer = buffer[positions[-1]:]
                if len(buffer) > 16*1024*1024:
                    raise UnsupportedFrontierLayout('individual record/layout exceeds buffer budget')
            offset += len(block)
        if buffer.strip():
            if index in wanted:
                found[index] = json.loads(buffer.strip().rstrip(b',').strip())
            index += 1
        if set(found) != wanted:
            raise ValueError('saved frontier refers to missing record indices')
    if (path.stat().st_size, path.stat().st_mtime_ns) != signature:
        raise ValueError('saved GP result changed during reading')
    completion = path.parent.parent / (data['galaxy']+'_completion.json')
    if completion.is_file():
        expected = json.loads(completion.read_text()).get('records')
        if expected is not None and expected != index:
            raise ValueError('candidate count differs from saved completion metadata')
    old_ids = sorted(found)
    remap = {old: i for i, old in enumerate(old_ids)}
    data['records'] = [dict(found[old], original_candidate_index=old) for old in old_ids]
    for item in data['frontier']:
        item['original_record_index'] = item['record_index']
        item['record_index'] = remap[item['record_index']]
    data['selection_reader'] = {'mode': 'frontier-stream', 'evaluated_records': index,
        'loaded_records': len(old_ids), 'frontier_logp_epsilon': 1e-12,
        'discarded_record_syntax_fully_validated': False,
        'saved_result_sha256': digest.hexdigest()}
    return data, digest.hexdigest()
