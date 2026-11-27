
from .matcher import scan_record

def scan_with_scheduler(records, guide, max_mismatches=3, use_dask=True):
    if use_dask:
        try:
            from dask import compute, delayed

            jobs = [
                delayed(scan_record)(r, guide, max_mismatches)
                for r in records
            ]

            groups = compute(*jobs, scheduler="threads")
            hits = [h for group in groups for h in group]

            return sorted(
                hits,
                key=lambda h: (h.mismatches, h.record_id, h.start)
            )

        except ImportError:
            pass

    hits = [
        h for r in records
        for h in scan_record(r, guide, max_mismatches)
    ]

    return sorted(
        hits,
        key=lambda h: (h.mismatches, h.record_id, h.start)
    )
