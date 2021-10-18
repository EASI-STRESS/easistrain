import configparser


def load_ini(filename):
    config = configparser.ConfigParser()
    config.read(filename)
    config = dict(config["arguments"])
    result = dict()
    for name, string_value in config.items():
        if string_value == "None":
            result[name] = None
        elif string_value.isdigit():
            result[name] = int(string_value)
        else:
            result[name] = string_value
    return result


if __name__ == "__main__":
    from pprint import pprint

    pprint(load_ini("exe_integration_2D.ini"))
