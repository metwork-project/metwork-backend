from enum import IntEnum


class AnnotationStatus(IntEnum):
    UNDEFINED = 0
    EXPLORED = 10
    PUTATIVE = 20
    VALIDATED = 30
    ERROR = 90
