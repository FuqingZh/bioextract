
from importlib import import_module
from types import ModuleType
from typing import TYPE_CHECKING, Any

__all__ = ["stringdb"]

if TYPE_CHECKING:
    import bioextract.stringdb as stringdb

_ALIAS_MODULES: dict[str, str] = {
    "stringdb": "bioextract.stringdb",
}


def __getattr__(name: str) -> Any:
    module_name = _ALIAS_MODULES.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_loaded: ModuleType = import_module(module_name)
    globals()[name] = module_loaded
    return module_loaded


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
