"""Self-registration registry for `ModelBackend` implementations.

To add a new experiment script as a backend: implement `ModelBackend` in a
new module under this package, then import that module from
`backend/model_backends/__init__.py` so it registers itself. Nothing else
needs to change.
"""

from __future__ import annotations

from backend.model_backends.base import ModelBackend

DEFAULT_BACKEND = "cross_attention"

_BACKENDS: dict[str, ModelBackend] = {}


def register(backend: ModelBackend) -> ModelBackend:
    _BACKENDS[backend.name] = backend
    return backend


def get_backend(name: str) -> ModelBackend:
    try:
        return _BACKENDS[name]
    except KeyError as exc:
        raise KeyError(f"Unknown model backend: {name!r}. Registered: {sorted(_BACKENDS)}") from exc
