from __future__ import annotations
from typing import Optional
from datetime import date
from matplotlib import pyplot as plt
from dataclasses import dataclass

PLOT_HIGHT = 5  # Высота графика
PLOT_WIDTH = 10  # Ширина графика


@dataclass
class PlotData:
    data: pd.Series
    style: str = ''
    width: int = 1
    label: str = ''


@dataclass
class Scale:
    min: int = 0
    max: int = 1


def build_plot(
    data: dict,
):
    """Построитель графиков"""

    _, ax = plt.subplots(figsize=(PLOT_WIDTH, PLOT_HIGHT))

    # Проинициализировать дополнительные шкалы
    axes = [ax]
    axes.extend([ax.twinx() for _ in range(len(data[1:]))])

    for _ax in axes[2:]:
        _ax.spines.right.set_position(("axes", 1.15))

    plots = []
    xlabel = '$Время, сутки$'

    # Построить графики
    for ax_index in range(len(data)):
        maxes = []
        mins = []

        _ax = axes[ax_index]

        x = data[ax_index].get('x')
        dt = data[ax_index].get('dt')
        mu = data[ax_index].get('mu')
        x_scale = data[ax_index].get('x_scale')
        y_scale = data[ax_index].get('y_scale')

        x_min = 0
        x_max = len(x)
        if x_scale:
            x_min = int(x_scale.min / dt)
            x_max = int(x_scale.max / dt)
        x = x[x_min:x_max]

        for plot_data in data[ax_index].get('data'):

            values = plot_data.data[x_min:x_max]

            (p,) = _ax.plot(x, values, plot_data.style, lw=plot_data.width, label=f"${plot_data.label}$")
            plots.append(p)

            maxes.append(values.max())
            mins.append(values.min())

        # Масштаб по измерениям
        min_limit = y_scale.min if y_scale else min(mins) * 0.9
        max_limit = y_scale.max if y_scale else max(maxes) * 1.1

        y_label = f"${', '.join([plot_data.label for plot_data in data[ax_index].get('data') if plot_data.label])}$"
        if mu:
            y_label += f', ${mu}$'

        _ax.set(ylim=(min_limit, max_limit), xlabel=xlabel, ylabel=y_label)

    ax.legend()
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
    #          fancybox=True, shadow=True, ncol=5, handles=plots)
    ax.grid()
    return plt


def get_dict_for_data():
    return {
        'q_N': [],
        'p_3': [],
        'p_1': [],
        'p_2': [],
        'q': [],
        'h_4': [],
        'u': [],
        'p_4': [],
        'p_8': [],
        'p_5': [],
        'N_1': [],
        'N_2': [],
        'n_NN': [],  # КПД ЭЦН
        'betta_G3': [],
        'betta_GN': [],
        'q_L': [],
        'agzu': [],  # Моменты работы АГЗУ
        'p_L': [],
        'b_betta_3': [],
        'b_betta_L': [],
        'x': [],
    }
