import os
from datetime import datetime


def log_parameters(filename, parameters, task_name):
    filename = os.path.join(parameters["root_data"], filename)
    with open(filename, "w") as fwlog:
        fwlog.write(f"{task_name.upper()} LOG FILE\n")
        fwlog.write(f"Date and time : {datetime.now()}\n")
        fwlog.write(
            f"#$#$#$#$#$#$#The arguments used for {task_name.lower()} are below: \n"
        )
        for name, value in parameters.items():
            fwlog.write(f"{name} = {value}\n")
        fwlog.write("************____________________________________**************\n")
