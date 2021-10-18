import json
from easistrain.load_ini import load_ini

with open("parameters.json", "w") as f:
    config = load_ini("exe_integration_2D.ini")
    json.dump(dict(config), f, indent=2)
