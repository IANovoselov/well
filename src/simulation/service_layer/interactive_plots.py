"""Интерактивные графики на Plotly.

Аналог build_plot() из services.py: принимает тот же формат plot_data.

    from src.simulation.service_layer.interactive_plots import build_plot_interactive
    build_plot_interactive(plot_data)

Зависимость: plotly (Jupyter).
"""
from __future__ import annotations

from typing import Any

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.simulation.service_layer.services import PLOT_HIGHT, PLOT_WIDTH, PlotData, Scale

X_LABEL = '$Время, сутки$'
PLOT_WIDTH_PX = int(PLOT_WIDTH * 96)
PLOT_HEIGHT_PX_PER_PANEL = int(PLOT_HIGHT * 96)
PLOT_BG = '#ffffff'
PAPER_BG = '#ffffff'
GRID_COLOR = '#cccccc'

LATEX_FONT: dict[str, Any] = {
    'family': 'STIXGeneral, "Computer Modern", "Latin Modern Math", serif',
    'size': 12,
    'color': '#000000',
}

_LINE_STYLES: dict[str, dict[str, Any]] = {
    '': {},
    'r': {'color': '#d62728'},
    'k': {'color': '#000000'},
    'g': {'color': '#2ca02c'},
    'k--': {'color': '#000000', 'dash': 'dash'},
    'g--': {'color': '#2ca02c', 'dash': 'dash'},
    'r--': {'color': '#d62728', 'dash': 'dash'},
}


def _as_array(series_or_list) -> pd.Series:
    if isinstance(series_or_list, pd.Series):
        return series_or_list
    return pd.Series(series_or_list)


def _time_to_index(time_days: float, dt: float) -> int:
    return max(0, int(time_days / dt))


def _line_style(plot_data: PlotData) -> dict[str, Any]:
    style = _LINE_STYLES.get(plot_data.style, {})
    return {'width': plot_data.width, **style}


def _math_text(text: str) -> str:
    if not text:
        return text
    stripped = text.strip()
    if stripped.startswith('$') and stripped.endswith('$'):
        return stripped
    return f'${stripped}$'


def _panel_ylabel(panel: dict) -> str:
    labels = [plot_data.label for plot_data in panel.get('data', []) if plot_data.label]
    y_label = f"${', '.join(labels)}$" if labels else ''
    mu = panel.get('mu')
    if mu:
        y_label += f', ${mu}$'
    return y_label


def _panel_slice_indices(panel: dict) -> tuple[int, int]:
    """Индексы окна по времени — как x_scale в build_plot()."""
    x = _as_array(panel['x'])
    dt = panel['dt']
    n_points = len(x)
    x_min_idx = 0
    x_max_idx = n_points
    x_scale: Scale | None = panel.get('x_scale')
    if x_scale:
        x_min_idx = _time_to_index(float(x_scale.min), dt)
        x_max_idx = _time_to_index(float(x_scale.max), dt)
    if x_max_idx <= x_min_idx:
        x_max_idx = min(n_points, x_min_idx + 1)
    return x_min_idx, min(n_points, x_max_idx)


def _slice_panel(panel: dict, x_min_idx: int, x_max_idx: int) -> tuple[pd.Series, list[tuple[PlotData, pd.Series]]]:
    x = _as_array(panel['x'])
    x_slice = x.iloc[x_min_idx:x_max_idx]
    traces = []
    for plot_data in panel.get('data', []):
        values = _as_array(plot_data.data).iloc[x_min_idx:x_max_idx]
        traces.append((plot_data, values))
    return x_slice, traces


def build_plotly_figure(data: list[dict]) -> go.Figure:
    """Построить Plotly-фигуру для plot_data."""
    if not data:
        raise ValueError('plot_data must contain at least one panel')

    fig = make_subplots(
        rows=len(data),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.06,
    )

    for row, panel in enumerate(data, start=1):
        x_min_idx, x_max_idx = _panel_slice_indices(panel)
        x_slice, traces = _slice_panel(panel, x_min_idx, x_max_idx)
        for plot_data, values in traces:
            fig.add_trace(
                go.Scatter(
                    x=x_slice,
                    y=values,
                    name=_math_text(plot_data.label) if plot_data.label else None,
                    line=_line_style(plot_data),
                    mode='lines',
                    hovertemplate=(
                        f'{_math_text(plot_data.label or "value")}'
                        '=%{y:.4f}<br>'
                        f'{_math_text("t")}=%{{x:.4f}} сут<extra></extra>'
                    ),
                ),
                row=row,
                col=1,
            )

        y_scale: Scale | None = panel.get('y_scale')
        if y_scale:
            fig.update_yaxes(range=[y_scale.min, y_scale.max], row=row, col=1)
        else:
            mins, maxes = [], []
            for _, values in traces:
                mins.append(float(values.min()))
                maxes.append(float(values.max()))
            if mins and maxes:
                fig.update_yaxes(
                    range=[min(mins) * 0.9, max(maxes) * 1.1],
                    row=row,
                    col=1,
                )

        fig.update_yaxes(
            title_text=_panel_ylabel(panel),
            title_standoff=8,
            title_font=LATEX_FONT,
            tickfont=LATEX_FONT,
            showgrid=True,
            gridcolor=GRID_COLOR,
            gridwidth=1,
            zeroline=False,
            row=row,
            col=1,
        )
        fig.update_xaxes(
            title_text=X_LABEL if row == len(data) else None,
            title_font=LATEX_FONT,
            tickfont=LATEX_FONT,
            showgrid=True,
            gridcolor=GRID_COLOR,
            gridwidth=1,
            zeroline=False,
            row=row,
            col=1,
        )

    fig.update_layout(
        width=PLOT_WIDTH_PX,
        height=PLOT_HEIGHT_PX_PER_PANEL * len(data),
        plot_bgcolor=PLOT_BG,
        paper_bgcolor=PAPER_BG,
        font=LATEX_FONT,
        hovermode='x unified',
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.02,
            xanchor='right',
            x=1,
            font=LATEX_FONT,
        ),
        margin=dict(t=60, b=55, l=90, r=40),
    )
    return fig


def build_plot_interactive(data: list[dict]) -> go.Figure:
    """Интерактивный аналог build_plot() на Plotly (зум, pan, hover)."""
    fig = build_plotly_figure(data)

    return fig
