# Идентификация через фильтр

import scipy
import pandas as pd
import numpy as np

from src.helpers import filter_data


def ident_values(ident_k, ident_dt, df_ident, well, filter_name=None):
  calc_q_L = []
  calc_dp = []
  calc_dp_t = []
  calc_q_t = []
  calc_p_1_t = []
  calc__q_t = []

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

    p_1_t = _p_3 + (well.oil.gamma/b_betta_3_t)*(well.params.H_R-well.pump.H_N+well.params.r_K*q_t)

    calc_q_L.append(q_L)
    calc_dp.append(dp)
    calc_dp_t.append(dp_t)
    calc_q_t.append(q_t)
    calc_p_1_t.append(p_1_t)
    calc__q_t.append(_q_t)

    dp_t = dp_t + (ident_dt/well.pump.t_N)*(dp - dp_t)
    _q_t = _q_t + (ident_dt/well.reservoir.T_2)*(q_t - _q_t)
    _p_3 = _p_3 + (ident_dt/well.pump.t_N)*(df_ident['p_3'][_k] - _p_3)
    b_betta_3_t = b_betta_3_t + (ident_dt/well.pump.t_N)*(b_betta_3 - b_betta_3_t)

  calc_df = pd.DataFrame({'q_L': calc_q_L,
                          'dp': calc_dp,
                          'dp_t': calc_dp_t,
                          'q_t': filter_data(calc_q_t, filter_name, ident_dt, well.pump.t_N),
                          'p_1_t': filter_data(calc_p_1_t, filter_name, ident_dt, well.pump.t_N),
                          '_q_t': filter_data(calc__q_t, filter_name, ident_dt, well.pump.t_N),
                          'q_t_original': calc_q_t,
                          'p_1_t_original': calc_p_1_t,
                          '_q_t_original': calc__q_t,
                          })

  return calc_df

def ident_values_nkt(ident_k, ident_dt, df_ident, well, filter_name=None):
  calc_q_N = []
  calc_b_betta_L = []
  calc_b_betta_3 = []
  calc_y = []
  calc_x_1 = []
  calc_x_2 = []
  calc_x_3 = []

  well.pump.i = 1

  for _k in range(ident_k):

      b_betta_L = well.b_0 - well.b_1*(well.params.p_G0 - df_ident['p_L'][_k])

      q_N = (b_betta_L/(well.oil.gamma*well.params.r_U))*(max(df_ident['p_8'][_k]-df_ident['p_L'][_k], 0))

      b_betta_3 = well.b_0 - well.params.alpha_G_3 * well.b_1 * (well.params.p_G0 - df_ident['p_3'][_k])

      calc_q_N.append(q_N)
      calc_b_betta_L.append(b_betta_L)
      calc_b_betta_3.append(b_betta_3)

      if _k == 0:
          _p_3 = df_ident['p_3'][0]
          _b_betta_L = (1/b_betta_L)
      else:
          _p_3 = _p_3 + (ident_dt / well.pump.t_N) * (df_ident['p_3'][_k] - _p_3)
          _b_betta_L = _b_betta_L + (ident_dt / well.pump.t_N) * ((1/b_betta_L) - _b_betta_L)

      y = _b_betta_L*well.pump.H_N + ((df_ident['p_8'][_k] - _p_3)/well.oil.gamma)
      calc_y.append(y)

      well.pump.update_segment(q_N)
      well.pump.update_segment(q_N)

      x_1 = (q_N/b_betta_L)
      x_2 = (well.pump.lambda_0[well.pump.i]*well.pump.h_0*df_ident['u'][_k] *df_ident['u'][_k] )/b_betta_3
      x_3 = ((well.pump.lambda_1[well.pump.i]*well.pump.h_0*df_ident['u'][_k] *q_N)/(well.pump.q_0 * b_betta_3))

      if _k == 0:
          _x_1 = x_1
          _x_2 = x_2
          _x_3 = x_3
      else:
          _x_1 = _x_1 + (ident_dt / well.pump.t_N) * (x_1 - _x_1)
          _x_2 = _x_2 + (ident_dt / well.pump.t_N) * (x_2 - _x_2)
          _x_3 = _x_3 + (ident_dt / well.pump.t_N) * (x_3 - _x_3)


      calc_x_1.append(_x_1)
      calc_x_2.append(_x_2)
      calc_x_3.append(_x_3)

  calc_df = pd.DataFrame({'q_N': calc_q_N,
                          'b_betta_L': calc_b_betta_L,
                          'b_betta_3': calc_b_betta_3,
                          'y': calc_y,
                          'x_1': calc_x_1,
                          'x_2': calc_x_2,
                          'x_3': calc_x_3,
                          })

  return calc_df

def ident_values_inflow(ident_k, ident_dt, df_ident, well, filter_name=None):
    calc_q_N = []
    calc_b_betta_L = []
    calc_b_betta_3 = []
    calc_y = []
    calc_x_1 = []
    calc_x_2 = []
    calc_x_3 = []

    t_t = well.pump.t_N * 5

    for _k in range(ident_k):

        b_betta_L = well.b_0 - well.b_1*(well.params.p_G0 - df_ident['p_L'][_k])

        q_N = (b_betta_L/(well.oil.gamma*well.params.r_U))*(df_ident['p_8'][_k]-df_ident['p_L'][_k])

        b_betta_3 = well.b_0 - well.params.alpha_G_3 * well.b_1 * (well.params.p_G0 - df_ident['p_3'][_k])

        d_p = df_ident['p_3'][_k]-df_ident['p_4'][_k]

        calc_q_N.append(q_N)
        calc_b_betta_L.append(b_betta_L)
        calc_b_betta_3.append(b_betta_3)

        if _k == 0:
            _q_N = q_N
            _p_3 = df_ident['p_3'][0]
            _b_betta_3 = (1/b_betta_3)
            _d_p = d_p
        else:
            _p_3 = _p_3 + (ident_dt / t_t) * (df_ident['p_3'][_k] - _p_3)
            _b_betta_3 = _b_betta_3 + (ident_dt / t_t) * ((1/b_betta_3) - _b_betta_3)
            _q_N = _q_N + (ident_dt / t_t) * (q_N - _q_N)
            _d_p = _d_p + (ident_dt / t_t) * (d_p - _d_p)

        q_tilda = _q_N + (well.params.S_t / (well.oil.gamma * t_t)) * (d_p - _d_p)

        if _k == 0:
            q_tilda_T = q_tilda
        else:
            q_tilda_T = q_tilda_T + (ident_dt / well.reservoir.T_2) * (q_tilda - q_tilda_T)

        y = _p_3 + well.oil.gamma * _b_betta_3 * (well.params.H_R - well.pump.H_N + well.params.r_K * q_tilda)
        calc_y.append(y)

        x_1 = 1
        x_2 = q_tilda
        x_3 = q_tilda_T

        calc_x_1.append(x_1)
        calc_x_2.append(x_2)
        calc_x_3.append(x_3)

    calc_df = pd.DataFrame({'q_N': calc_q_N,
                            'b_betta_L': calc_b_betta_L,
                            'b_betta_3': calc_b_betta_3,
                            'y': calc_y,
                            'x_1': calc_x_1,
                            'x_2': calc_x_2,
                            'x_3': calc_x_3,
                            })

    return calc_df

def identificate_nkt(calc_df, ident_k):
  X = np.array([-calc_df['x_1'],
                calc_df['x_2'],
                -calc_df['x_3']]).T
  y = np.array(calc_df['y']).T
  b, squared_error_sum, matrix_rank, SVD_ = scipy.linalg.lstsq(X, y)

  #b = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y)

  return b, squared_error_sum

def identificate_inflow(calc_df, ident_k):
  X = np.array([calc_df['x_1'],
                -calc_df['x_2'],
                -calc_df['x_3']]).T
  y = np.array(calc_df['y']).T
  b, squared_error_sum, matrix_rank, SVD_ = scipy.linalg.lstsq(X, y)

  #b = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y)

  return b, squared_error_sum

def identificate(calc_df, ident_k):
  X = np.array([[1]*ident_k,
              calc_df['q_t']*-1,
              calc_df['_q_t']*-1]).T
  y = np.array(calc_df['p_1_t']).T
  b, squared_error_sum, matrix_rank, SVD_ = scipy.linalg.lstsq(X, y)

  #b = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y)

  return b, squared_error_sum

def get_data_by_slices(calc_df, ident_dt, times):

  slices = []

  for start, end in times:
    start_index = int(start / ident_dt)
    end_index = int(end / ident_dt)

    slices.append(calc_df.iloc[start_index:end_index])

  result = pd.concat(slices, ignore_index=True)

  return result

def identificate_regul(calc_df, calc_df_static, v_r, m2=0, m3=0):
    # Динамический участок
    X1 = np.column_stack([
        np.ones(len(calc_df)),
        -calc_df['q_t'].to_numpy(),
        -calc_df['_q_t'].to_numpy(),
    ])
    y1 = calc_df['p_1_t'].to_numpy()
    # Статический участок
    X2 = np.column_stack([
        np.ones(len(calc_df_static)),
        -calc_df_static['_q_t'].to_numpy(),
        -calc_df_static['_q_t'].to_numpy(),
    ])
    y2 = calc_df_static['p_1_t'].to_numpy()


    blocks_X = [X1, m2 * X2]
    blocks_y = [y1, m2 * y2]
    # Регуляризация: r_1 ≈ v_r * r_2
    if m3 > 0:
        X3 = np.sqrt(m3) * np.array([[0.0, 1.0, -v_r]])
        y3 = np.array([0.0])
        blocks_X.append(X3)
        blocks_y.append(y3)
    X = np.vstack(blocks_X)
    y = np.concatenate(blocks_y)
    b, squared_error_sum, matrix_rank, SVD_ = scipy.linalg.lstsq(X, y)
    return b

    # формируем и заполняем матрицу размерностью 2x2
    A1 = np.empty((3, 3))
    calc_df_len = len(calc_df)
    A1[[0], [0]] = calc_df_len
    A1[[0], [1]] = -calc_df['q_t'].sum()
    A1[[0], [2]] = -calc_df['_q_t'].sum()

    A1[[1], [0]] = -calc_df['q_t'].sum()
    A1[[1], [1]] = sum(value**2 for value in calc_df['q_t'])
    A1[[1], [2]] = sum(calc_df['_q_t'][i]*calc_df['q_t'][i] for i in range(calc_df_len))

    A1[[2], [0]] = -calc_df['_q_t'].sum()
    A1[[2], [1]] = sum(calc_df['_q_t'][i]*calc_df['q_t'][i] for i in range(calc_df_len))
    A1[[2], [2]] = sum(value**2 for value in calc_df['_q_t'])

    m1 = 1
    A1 = m1 * A1

    calc_df_static_len = len(calc_df_static)
    A2 = np.empty((3, 3))
    A2[[0], [0]] = calc_df_static_len
    A2[[0], [1]] = -calc_df_static['_q_t'].sum()
    A2[[0], [2]] = -calc_df_static['_q_t'].sum()

    A2[[1], [0]] = -calc_df_static['_q_t'].sum()
    A2[[1], [1]] = sum([value**2 for value in calc_df_static['_q_t']])
    A2[[1], [2]] = sum([value**2 for value in calc_df_static['_q_t']])

    A2[[2], [0]] = -calc_df_static['_q_t'].sum()
    A2[[2], [1]] = sum([value**2 for value in calc_df_static['_q_t']])
    A2[[2], [2]] = sum([value**2 for value in calc_df_static['_q_t']])

    m2 = 1
    A2 = m2 * A2

    A3 = np.empty((3, 3))
    A3[[0], [0]] = 0
    A3[[0], [1]] = 0
    A3[[0], [2]] = 0

    A3[[1], [0]] = 0
    A3[[1], [1]] = 1
    A3[[1], [2]] = -v_r

    A3[[2], [0]] = 0
    A3[[2], [1]] = -v_r
    A3[[2], [2]] = v_r * v_r

    A3 = m3 * A3

    A = A1 + A2 + A3

    # находим обратную матрицу
    A = np.linalg.inv(A)
    # формируем и заполняем матрицу размерностью 3x1
    C1 = np.empty((3, 1))
    C1[0] = calc_df['p_1_t'].sum()
    C1[1] = -sum(calc_df['p_1_t'][i]*calc_df['q_t'][i] for i in range(calc_df_len))
    C1[2] = -sum(calc_df['p_1_t'][i]*calc_df['_q_t'][i] for i in range(calc_df_len))

    C2 = np.empty((3, 1))
    C2[0] = calc_df_static['p_1_t'].sum()
    C2[1] = -sum([calc_df_static['p_1_t'][i]*calc_df_static['_q_t'][i] for i in range(calc_df_static_len)])
    C2[2] = -sum([calc_df_static['p_1_t'][i]*calc_df_static['_q_t'][i] for i in range(calc_df_static_len)])

    C2 = m2 * C2

    C = C1 + C2

    # умножаем матрицу на вектор
    ww = np.dot(A, C)
    return ww


def get_res_dict_for_ident(quant_num, denominator_num, noise_num, experiment_num):
    return {'p_R': np.zeros((quant_num, denominator_num, noise_num, experiment_num)),
            'r_1': np.zeros((quant_num, denominator_num, noise_num, experiment_num)),
            'r_2': np.zeros((quant_num, denominator_num, noise_num, experiment_num)),
             }
