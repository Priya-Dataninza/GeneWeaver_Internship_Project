
try:
    from numba import cuda
    NUMBA_AVAILABLE = True
except Exception:
    cuda = None
    NUMBA_AVAILABLE = False

def cuda_is_available():
    if not NUMBA_AVAILABLE:
        return False

    try:
        return bool(cuda.is_available())
    except Exception:
        return False

def choose_engine(requested="auto"):
    if requested == "cpu":
        return "cpu"

    if requested == "gpu":
        if not cuda_is_available():
            raise RuntimeError("GPU requested but CUDA unavailable")
        return "gpu"

    if requested == "auto":
        return "gpu" if cuda_is_available() else "cpu"

    raise ValueError("Use auto, cpu or gpu")