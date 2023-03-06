from ewokscore import load_graph


def create_workflow(filename):
    # Define a workflow
    nodes = [
        {
            "id": "Integrate2D",
            "class": "easistrain.task_integration_2D.Integrate2D",
        },
    ]
    # links = [
    #    {"source": "task1", "target": "task2", "arguments": {"a": "result"}},
    #    {"source": "task2", "target": "task3", "arguments": {"a": "result"}},
    # ]
    workflow = {"nodes": nodes, "links": []}

    # Execute a workflow (use a proper Ewoks task scheduler in production)
    graph = load_graph(workflow)
    graph.dump(filename, indent=2)
