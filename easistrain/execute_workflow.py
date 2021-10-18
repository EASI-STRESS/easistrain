from ewokscore import load_graph
import json
from easistrain.log_parameters import log_parameters


def execute_graph1(workflow_filename, param_filename):
    # Load parameters
    with open(param_filename, "r") as f:
        parameters = json.load(f)

    log_parameters("exe_integration.log", parameters, "2d integration")

    # Load graph
    graph = load_graph(workflow_filename)

    # Loop inside Integrate2D
    graph.graph.nodes["Integrate2D"]["inputs"] = parameters

    varinfo = {"root_uri": None}  # optional
    tasks = graph.execute(varinfo=varinfo)
    for name, task in tasks.items():
        print(name, task.output_values)


def execute_graph2(workflow_filename, param_filename):
    # Load parameters
    with open(param_filename, "r") as f:
        parameters = json.load(f)

    log_parameters("exe_integration.log", parameters, "2d integration")

    # Load graph
    graph = load_graph(workflow_filename)

    # Loop outside Integrate2D
    numscanstart, numscanend = parameters.pop("numScan")

    # Execute graph for each scan
    for numscan in range(numscanstart, numscanend + 1):
        parameters["numScan"] = [numscan, numscan + 1]
        graph.graph.nodes["Integrate2D"]["inputs"] = parameters

        varinfo = {"root_uri": None}  # optional
        tasks = graph.execute(varinfo=varinfo)
        for name, task in tasks.items():
            print(name, task.output_values)
