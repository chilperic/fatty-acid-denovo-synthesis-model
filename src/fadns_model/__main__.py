"""Run an explicit JSON configuration and export results with provenance."""

import argparse
import json
from pathlib import Path

import numpy as np

from .provenance import write_manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if args.output.exists():
        parser.error("Output path already exists; choose a new run directory")
    from . import basic
    from .simulation import simulate_time_course
    from .variants import simulate_variant

    variant = config["variant"]
    if variant == "basic_model":
        basic.validate_parameters(config["parameters"])
        result = simulate_time_course(
            basic.rhs, config["initial_state"], config["parameters"], **config["solver"]
        )
        compounds = basic.COMPOUNDS
    else:
        result = simulate_variant(
            variant, config["initial_state"], config["parameters"], **config["solver"]
        )
        compounds = result.compounds
    args.output.mkdir(parents=True)
    np.savetxt(
        args.output / "trajectory.csv",
        np.column_stack([result.t, result.y.T]),
        delimiter=",",
        header=",".join(["time", *compounds]),
        comments="",
    )
    (args.output / "solver.json").write_text(
        json.dumps(
            {
                "success": bool(result.success),
                "message": result.message,
                "nfev": int(result.nfev),
                "njev": int(result.njev),
            },
            indent=2,
        )
        + "\n"
    )
    write_manifest(args.output, config, args.config, "fadns-model", ["numpy", "scipy"])


if __name__ == "__main__":
    main()
