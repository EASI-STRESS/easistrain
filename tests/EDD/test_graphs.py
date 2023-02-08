from pathlib import Path

from ewoksutils.import_utils import qualname
from ewokscore import execute_graph

from easistrain.EDD.calibrationEDD import calibEdd
from .test_calibration import calib_edd_init, calib_edd_assert


def edd_graph(config):
    default_inputs = [{"name": name, "value": value} for name, value in config.items()]
    nodes = [
        {
            "id": "calib",
            "task_type": "method",
            "task_identifier": qualname(calibEdd),
            "default_inputs": default_inputs,
        }
    ]
    links = []
    return {"nodes": nodes, "links": links}


def test_edd_graph(tmp_path: Path):
    test_data_path, config = calib_edd_init(tmp_path)
    graph = edd_graph(config)
    results = execute_graph(graph)
    assert results == {"return_value": None}

    calib_edd_assert(test_data_path, config)
