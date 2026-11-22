from dataclasses import dataclass


@dataclass
class MatchResult:
    record_id: str
    start: int
    matched_sequence: str
    mismatches: int
    mismatch_positions: tuple
