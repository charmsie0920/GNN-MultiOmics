"""Backend wrapper for experiments/Cross Attention Fusion/cross_attention_baseline.py.

Runs the script unmodified by loading it as a module and monkeypatching its
`print` (to also stream lines to the UI) and its `train` function (to stop
streaming once training finishes, so the final results table the script
prints afterwards only reaches the real terminal, not the UI console).
"""

from __future__ import annotations

import builtins
import importlib.util
import types
from collections.abc import Callable
from pathlib import Path

from backend.model_backends.base import ModelBackend
from backend.model_backends.registry import register

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_PATH = _REPO_ROOT / "experiments" / "Cross Attention Fusion" / "cross_attention_baseline.py"


def _load_module() -> types.ModuleType:
    spec = importlib.util.spec_from_file_location("cross_attention_baseline_run", _SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module spec from {_SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CrossAttentionBackend(ModelBackend):
    name = "cross_attention"
    # Calibrated from the user's measured runs: ~7 min on CUDA, ~30 min on CPU.
    expected_duration_seconds = {"cuda": 420.0, "cpu": 1800.0}

    def run(self, log: Callable[[str], None]) -> None:
        module = _load_module()
        streaming = {"enabled": True}

        def tee_print(*args: object, sep: str = " ", end: str = "\n", **kwargs: object) -> None:
            builtins.print(*args, sep=sep, end=end, **kwargs)
            if streaming["enabled"]:
                log(sep.join(str(arg) for arg in args))

        module.print = tee_print

        original_train = module.train

        def wrapped_train(*args: object, **kwargs: object):
            result = original_train(*args, **kwargs)
            streaming["enabled"] = False
            return result

        module.train = wrapped_train

        module.main()


register(CrossAttentionBackend())
