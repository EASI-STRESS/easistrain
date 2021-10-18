import os
from easistrain.create_workflow import create_workflow
from easistrain.execute_workflow import execute_graph1


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Create and execute graphs")

    parser.add_argument(
        "--workflow",
        type=str,
        default="workflow.json",
        help="Workflow to execute",
    )

    parser.add_argument(
        "--parameters",
        type=str,
        default="parameters.json",
        help="Parameters for the workflow",
    )

    args, unknown = parser.parse_known_args()
    if not os.path.isfile(args.workflow):
        create_workflow(args.workflow)
    if not os.path.isfile(args.parameters):
        raise RuntimeError(f"Parameters file '{args.parameters}' does not exist")
    execute_graph1(args.workflow, args.parameters)
