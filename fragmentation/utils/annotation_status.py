from enum import IntEnum


class AnnotationStatus(IntEnum):
    UNDEFINED = 0
    PUTATIVE = 10
    VALIDATED = 20
    ERROR = 90

