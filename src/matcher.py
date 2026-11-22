from .models import MatchResult


def scan_record(record, guide, max_mismatches=3):
    guide = guide.upper()
    size = len(guide)
    hits = []

    for start in range(max(0, len(record.sequence) - size + 1)):
        window = record.sequence[start:start + size]

        positions = tuple(
            i for i, (a, b) in enumerate(zip(guide, window))
            if a != b
        )

        if len(positions) <= max_mismatches:
            hits.append(
                MatchResult(
                    record.record_id,
                    start,
                    window,
                    len(positions),
                    positions
                )
            )

    return hits


def scan_records(records, guide, max_mismatches=3):
    hits = [
        h
        for r in records
        for h in scan_record(r, guide, max_mismatches)
    ]

    return sorted(
        hits,
        key=lambda h: (h.mismatches, h.record_id, h.start)
    )
