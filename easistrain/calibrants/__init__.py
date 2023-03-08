import os

CALIBRANT_ROOT_DIR = os.path.abspath(os.path.dirname(__file__))


def calibrant_filename(filename: str) -> str:
    if os.path.exists(filename):
        return filename

    local_filename = os.path.join(CALIBRANT_ROOT_DIR, os.path.basename(filename))
    if os.path.exists(local_filename):
        return local_filename

    if not os.path.splitext(local_filename)[-1]:
        local_filename += ".dat"

    if os.path.exists(local_filename):
        return local_filename

    raise FileNotFoundError(filename)
