from ewokscore import load_graph

# Define a workflow
nodes = [
    {
        "id": "Integrate2D",
        "class": "easistrrain.class_integrate_2D.Integrate2D",
        "inputs": {
            "root_data": "/home/esrf/slim/data/blc_13060/TiC_GG/",
            "h5file": "blc13060_TiC_GG.h5",
            "scan": "TiC_GG_0002_1.1",
            "detector_name": "p3",
            "poni_file": "/home/esrf/slim/data/blc_13060/CeO2/CeO2_image.poni",
            "npt_rad": "5000",
            "npt_azim": "72",
            "x_unit": "2th_deg",
            "im_dark": "0",
            "im_mask": "/home/esrf/slim/data/blc_13060/CeO2/CeO2_image_mask.edf",
        },
	{
        "id": "fit",
        "class": "easistrrain.class_fitting_peaks.Fit",
        "inputs": {
            "...
        },
    },
]
# links = [
#    {"source": "task1", "target": "task2", "arguments": {"a": "result"}},
#    {"source": "task2", "target": "task3", "arguments": {"a": "result"}},
# ]
workflow = {"nodes": nodes, "links": []}

# Execute a workflow (use a proper Ewoks task scheduler in production)
graph = load_graph(workflow)
graph.dump("test.json", indent=2)
