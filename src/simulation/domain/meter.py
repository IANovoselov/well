import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from datetime import datetime, timedelta
from scipy.stats import truncnorm
import pandas as pd


class Meter:

  def __init__(self, dt=1):
    self.noise_enable = False
    self.noise_amplitude = 0.01
    self.dt = dt

  def measure(self, df, quant_step=0.1, denominator=7):
    """Сделать выборку измерений"""

    ident_k = (len(df)//denominator)
    ident_dt = self.dt * denominator

    if quant_step is not None:

      df_ident = pd.DataFrame({'p_3': [np.round(index / quant_step) * quant_step for index in self.add_noise(df['p_3'])[::denominator]],
                               'p_4': [np.round(index / quant_step) * quant_step for index in self.add_noise(df['p_4'])[::denominator]],
                               'p_8': [np.round(index / quant_step) * quant_step for index in self.add_noise(df['p_8'])[::denominator]],
                               'p_L': [np.round(index / quant_step) * quant_step for index in self.add_noise(df['p_L'])[::denominator]],
                               'x': [_x for _x in df['x'][::denominator]],
                               })
    else:
      df_ident = pd.DataFrame({'p_3': [_x for _x in self.add_noise(df['p_3'])[::denominator]],
                               'p_4': [_x for _x in self.add_noise(df['p_4'])[::denominator]],
                               'p_8': [_x for _x in self.add_noise(df['p_8'])[::denominator]],
                               'p_L': [_x for _x in self.add_noise(df['p_L'])[::denominator]],
                               'x': [_x for _x in df['x'][::denominator]],
                               })

    return ident_dt, ident_k, df_ident

  def add_noise(self, data):
    if self.noise_enable and self.noise_amplitude:
      return add_noise(data, dt=self.dt, hf_amplitude=self.noise_amplitude)[0]

    return data


def add_noise(clean_data,
              dt,
              lf_amplitude=0.02,
              hf_amplitude=0.01,
              lf_cutoff=0.5):
    """
    Добавляет низкочастотный и высокочастотный шум к данным давления с привязкой к суточному времени

    Параметры:
    clean_data - массив чистых данных давления
    sample_interval_min - интервал между измерениями в минутах
    lf_amplitude - амплитуда низкочастотного шума (доля от диапазона данных)
    hf_amplitude - амплитуда высокочастотного шума (доля от диапазона данных)
    lf_cutoff - частота среза для НЧ шума
    random_seed - зерно для генератора случайных чисел

    Возвращает:
    noisy_data - зашумленные данные
    lf_noise - низкочастотный шум
    hf_noise - высокочастотный шум
    """

    #data_range = np.max(clean_data) - np.min(clean_data)
#
    #fs = 1.0 / dt  # частота дискретизации
#
    ## 1. Генерация низкочастотного шума (дрейф)
    ## Для суточных данных НЧ шум должен быть очень низкочастотным
    #b, a = signal.butter(2, lf_cutoff / (fs/2), 'low')
    #lf_noise = signal.lfilter(b, a, np.random.randn(n))
    #lf_noise = lf_amplitude * data_range * (lf_noise / np.max(np.abs(lf_noise)))

    # 2. Генерация высокочастотного шума
    hf_noise = generate_hf_noise(hf_amplitude*1.5, len(clean_data))

    # 3. Комбинация шумов с чистыми данными
    noisy_data = clean_data + hf_noise #+ lf_noise

    return noisy_data, hf_noise


def generate_hf_noise(amplitude: float, size: int) -> list[float]:
    """ Генерация высокочастотного шума """
    mu = 0
    sigma = amplitude
    min_val = -amplitude
    max_val = amplitude

    # Преобразование границ в параметры для усеченного распределения
    a = (min_val - mu) / sigma
    b = (max_val - mu) / sigma

    return truncnorm.rvs(a, b, loc=mu, scale=sigma, size=size)
