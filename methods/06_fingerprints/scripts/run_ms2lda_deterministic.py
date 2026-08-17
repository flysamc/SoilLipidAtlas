#!/usr/bin/env python3
"""Run the recovered MS2LDA producer with a fixed seed and one train worker."""

from __future__ import annotations

import json
import os
import runpy
import sys
from pathlib import Path

import MS2LDA.modeling as modeling


SEED = 20260804
TRAIN_WORKERS = 1


def output_directory(argv: list[str]) -> Path:
    if "--output" not in argv:
        return Path("ms2lda_results")
    position = argv.index("--output")
    if position + 1 >= len(argv):
        raise ValueError("--output requires a directory")
    return Path(argv[position + 1])


def main() -> None:
    producer = Path(os.environ["MS2LDA_BASE_PRODUCER"]).resolve()
    base_define_model = modeling.define_model
    base_train_model = modeling.train_model

    def deterministic_define_model(n_motifs, model_parameters=None):
        parameters = dict(model_parameters or {})
        parameters["seed"] = SEED
        return base_define_model(n_motifs, model_parameters=parameters)

    def deterministic_train_model(
        model,
        documents,
        iterations=100,
        train_parameters=None,
        convergence_parameters=None,
    ):
        parameters = dict(train_parameters or {})
        parameters["workers"] = TRAIN_WORKERS
        kwargs = {
            "iterations": iterations,
            "train_parameters": parameters,
        }
        if convergence_parameters is not None:
            kwargs["convergence_parameters"] = convergence_parameters
        return base_train_model(model, documents, **kwargs)

    modeling.define_model = deterministic_define_model
    modeling.train_model = deterministic_train_model
    runpy.run_path(str(producer), run_name="__main__")

    parameters_path = output_directory(sys.argv) / "run_parameters.json"
    parameters = json.loads(parameters_path.read_text(encoding="utf-8"))
    parameters["random_seed"] = SEED
    parameters["train_workers"] = TRAIN_WORKERS
    parameters["base_producer"] = str(producer)
    parameters["deterministic_wrapper"] = str(Path(__file__).resolve())
    parameters_path.write_text(json.dumps(parameters, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
