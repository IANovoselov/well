# Идентификация через фильтр

from dataclasses import dataclass, field
import scipy
import pandas as pd
import numpy as np

from scipy.signal import medfilt


@dataclass
class CalculatedValues:
  values: list[float] = field(default_factory=list)

  def add(self, element):
    self.values.append(element)

def ident_values(ident_k, ident_dt, df_ident, well, s_times=0):
  calc_q_L = CalculatedValues()
  calc_dp = CalculatedValues()
  calc_dp_t = CalculatedValues()
  calc_q_t = CalculatedValues()
  calc_p_1_t = CalculatedValues()
  calc__q_t = CalculatedValues()

  for _k in range(ident_k):

    b_betta_L = well.b_0 - well.b_1*(well.params.p_G0 - df_ident['p_L'][_k])

    q_L = (b_betta_L/(well.oil.gamma*well.params.r_U))*(df_ident['p_8'][_k]-df_ident['p_L'][_k])

    dp = (df_ident['p_3'][_k]-df_ident['p_4'][_k])

    b_betta_3 = well.b_0 - well.params.alpha_G_3 * well.b_1 * (well.params.p_G0 - df_ident['p_3'][_k])

    if _k == 0:
      dp_t = dp
      b_betta_3_t = b_betta_3

    q_t = q_L + ((well.params.S_t/(well.oil.gamma*well.pump.t_N)) * (dp - dp_t))

    if _k == 0:
      _q_t = 90
      _p_3 = df_ident['p_3'][0]

    if s_times and _k <= s_times / ident_dt:
        q_t = _q_t

    p_1_t = _p_3 + (well.oil.gamma/b_betta_3_t)*(well.params.H_R-well.pump.H_N+well.params.r_K*q_t)

    calc_q_L.add(q_L)
    calc_dp.add(dp)
    calc_dp_t.add(dp_t)
    calc_q_t.add(q_t)
    calc_p_1_t.add(p_1_t)
    calc__q_t.add(_q_t)

    dp_t = dp_t + (ident_dt/well.pump.t_N)*(dp - dp_t)
    _q_t = _q_t + (ident_dt/well.reservoir.T_2)*(q_t - _q_t)
    _p_3 = _p_3 + (ident_dt/well.pump.t_N)*(df_ident['p_3'][_k] - _p_3)
    b_betta_3_t = b_betta_3_t + (ident_dt/well.pump.t_N)*(b_betta_3 - b_betta_3_t)

  calc_df = pd.DataFrame({'q_L': calc_q_L.values,
                          'dp': calc_dp.values,
                          'dp_t': calc_dp_t.values,
                          'q_t': calc_q_t.values,
                          'p_1_t': calc_p_1_t.values,
                          '_q_t': calc__q_t.values,
                          })

  return calc_df

def identificate(calc_df, ident_k):
  X = np.array([[1 for _ in range(ident_k)],
              calc_df['q_t']*-1,
              calc_df['_q_t']*-1]).T
  y = np.array(calc_df['p_1_t']).T
  b, squared_error_sum, matrix_rank, SVD_ = scipy.linalg.lstsq(X, y)
  return b

def get_static_data(calc_df, ident_k, denominator, s_times, dt):

  s_q_t = []
  s_p_t = []

  for start, end in s_times:
    start_index = int((start / dt) // denominator)
    end_index = int((end / dt) // denominator)

    s_q_t.append(calc_df['_q_t'][start_index:end_index])
    s_p_t.append(calc_df['p_1_t'][start_index:end_index])

  s_q_t = pd.concat(s_q_t)
  s_p_t = pd.concat(s_p_t)

  k_s = len(s_q_t)

  m2 = ident_k / k_s
  return s_q_t, s_p_t, m2, k_s

def identificate_regul(calc_df, ident_k, denominator, v_r):

  # формируем и заполняем матрицу размерностью 2x2
  A1 = np.empty((3, 3))
  A1[[0], [0]] = ident_k
  A1[[0], [1]] = -sum(calc_df['q_t'])
  A1[[0], [2]] = -sum(calc_df['_q_t'])

  A1[[1], [0]] = -sum(calc_df['q_t'])
  A1[[1], [1]] = sum(value**2 for value in calc_df['q_t'])
  A1[[1], [2]] = sum(calc_df['_q_t'][i]*calc_df['q_t'][i] for i in range(ident_k))

  A1[[2], [0]] = -sum(calc_df['_q_t'])
  A1[[2], [1]] = sum(calc_df['_q_t'][i]*calc_df['q_t'][i] for i in range(ident_k))
  A1[[2], [2]] = sum(value**2 for value in calc_df['_q_t'])

  #A2 = np.empty((3, 3))
  #A2[[0], [0]] = k_s
  #A2[[0], [1]] = -sum(s_q_t)
  #A2[[0], [2]] = -sum(s_q_t)
#
  #A2[[1], [0]] = -sum(s_q_t)
  #A2[[1], [1]] = sum([value**2 for value in s_q_t])
  #A2[[1], [2]] = sum([value**2 for value in s_q_t])
#
  #A2[[2], [0]] = -sum(s_q_t)
  #A2[[2], [1]] = sum([value**2 for value in s_q_t])
  #A2[[2], [2]] = sum([value**2 for value in s_q_t])
#
  #A2 = m2 * A2

  A3 = np.empty((3, 3))
  A3[[0], [0]] = 0
  A3[[0], [1]] = 0
  A3[[0], [2]] = 0
#
  A3[[1], [0]] = 0
  A3[[1], [1]] = 1
  A3[[1], [2]] = -v_r
#
  A3[[2], [0]] = 0
  A3[[2], [1]] = -v_r
  A3[[2], [2]] = v_r * v_r

  A3 = ident_k * A3

  A = A1 + A3 #+ A2

  # находим обратную матрицу
  A = np.linalg.inv(A)
  # формируем и заполняем матрицу размерностью 3x1
  C1 = np.empty((3, 1))
  C1[0] = sum(calc_df['p_1_t'])
  C1[1] = -sum(calc_df['p_1_t'][i]*calc_df['q_t'][i] for i in range(ident_k))
  C1[2] = -sum(calc_df['p_1_t'][i]*calc_df['_q_t'][i] for i in range(ident_k))

  #C2 = np.empty((3, 1))
  #C2[0] = sum(s_p_t)
  #C2[1] = -sum([s_p_t[i]*s_q_t[i] for i in range(k_s)])
  #C2[2] = -sum([s_p_t[i]*s_q_t[i] for i in range(k_s)])
#
  #C2 = m2 * C2

  C = C1 #+ C2

  # умножаем матрицу на вектор
  ww = np.dot(A, C)
  return ww


def get_res_dict_for_ident(quant_num, denominator_num, noise_num, experiment_num):
    return {'p_R': np.zeros((quant_num, denominator_num, noise_num, experiment_num)),
            'r_1': np.zeros((quant_num, denominator_num, noise_num, experiment_num)),
            'r_2': np.zeros((quant_num, denominator_num, noise_num, experiment_num)),
             }