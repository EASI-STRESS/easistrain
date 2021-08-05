from ewokscore import load_graph

graph = load_graph("test.json")

#graph.nodes["Integrate2D"]["inputs"]["root_data"] = ...

varinfo = {"root_uri": None}  # optional
tasks = graph.execute(varinfo=varinfo)
for name,task in tasks.items():
    print(name, task.output_values)


