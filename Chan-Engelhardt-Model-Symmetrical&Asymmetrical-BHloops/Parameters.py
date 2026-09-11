
from pathlib import Path

# ============================================================
# READ PARAMETERS FROM FILE
# ============================================================

parameters = {}

file_path = Path(__file__).with_name("Parameters.txt")

with open(file_path, "r", encoding="utf-8") as f:
    for line in f:

        line = line.strip()

        # Skip empty lines and comments
        if not line or line.startswith("#"):
            continue

        name, value = line.split("=", 1)

        parameters[name.strip()] = float(value.strip())


# ============================================================
# GET PARAMETER BY NAME
# ============================================================

def get_parameter(name):
    if name not in parameters:
        raise KeyError(
            f"Parameter '{name}' is not found in Parameters.txt"
        )

    return parameters[name]