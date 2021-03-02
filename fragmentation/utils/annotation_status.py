from enum import IntEnum


class AnnotationStatus(IntEnum):
    UNDEFINED = 0
    UNRATED = 10
    PUTATIVE = 20
    VALIDATED = 30
    ERROR = 90


STATUS_LEVEL_MAPPING = {
    1: AnnotationStatus.VALIDATED,
    2: AnnotationStatus.PUTATIVE,
    3: AnnotationStatus.UNRATED,
}

LEVEL_STATUS_MAPPING = {v: k for k, v in STATUS_LEVEL_MAPPING.items()}
