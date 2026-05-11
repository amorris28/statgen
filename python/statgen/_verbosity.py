import sys

_LEVELS = ("quiet", "info")
_VERBOSITY = "info"


def set_verbosity(level: str) -> None:
    global _VERBOSITY
    if level not in _LEVELS:
        raise ValueError("verbosity must be one of: quiet, info")
    _VERBOSITY = level


def get_verbosity() -> str:
    return _VERBOSITY


def info(message: str) -> None:
    if _VERBOSITY == "info":
        print(message, file=sys.stderr)
