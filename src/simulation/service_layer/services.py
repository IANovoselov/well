from __future__ import annotations
from typing import Optional
from datetime import date

from simulation.domain.model import FirstOrderLag


def get_dict_for_data():
  return {'q_N': [],
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
                'x': []
                }