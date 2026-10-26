"""Public exception hierarchy for FastNSA."""


class FastNSAError(Exception):
    """Base class for errors raised by FastNSA."""


class FastNSAValidationError(FastNSAError, ValueError):
    """Raised when analysis inputs or parameters are invalid."""


class NativeBackendError(FastNSAError, RuntimeError):
    """Raised when the native FastNSA backend reports a failure."""


class FastNSACancelledError(NativeBackendError):
    """Raised when a native operation is canceled cooperatively."""


def backend_error(message: str) -> NativeBackendError:
    """Translate a native error message into a stable public exception."""
    if "cancel" in message.casefold():
        return FastNSACancelledError(message)
    return NativeBackendError(message)
