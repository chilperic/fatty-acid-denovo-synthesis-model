"""Plotting helpers."""

from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt


def plot_time_course(result, labels=None, output: str | Path | None = None):
    labels = labels or [f"state_{i}" for i in range(result.y.shape[0])]
    fig, ax = plt.subplots()
    for i, label in enumerate(labels):
        ax.plot(result.t, result.y[i], label=label)
    ax.set_xlabel("Time")
    ax.set_ylabel("Concentration / amount")
    ax.legend(frameon=False)
    fig.tight_layout()
    if output is not None:
        fig.savefig(output, bbox_inches="tight")
    return fig, ax
