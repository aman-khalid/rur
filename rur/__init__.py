from importlib import import_module

__all__ = ("uri", "uhmi", "painter", "utool", "drawer", "vr")


def __getattr__(name):
    if name in __all__:
        module = import_module(f"rur.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
