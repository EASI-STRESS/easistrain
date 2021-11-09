import yaml


def read_config_file(path: str):
    with open(path, "r") as config_file:
        return yaml.load(config_file, Loader=yaml.SafeLoader)
