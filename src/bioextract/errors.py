"""Stable exception categories shared by bioextract database handles."""


class BioextractError(RuntimeError):
    """Base class for errors raised through bioextract public APIs."""


class CapabilityError(BioextractError):
    """Raised when an input or publication lacks a requested capability."""


class IntegrityError(BioextractError):
    """Raised when source or publication integrity cannot be established."""


__all__ = [
    "BioextractError",
    "CapabilityError",
    "IntegrityError",
]
