import numpy as np

def calc_sco(p_R, r_1, r_2, reservoir):
  res = 0

  res += ((p_R/reservoir.p_R)-1)**2
  res += ((r_1/reservoir.r_1)-1)**2
  res += ((r_2/reservoir.r_2)-1)**2

  return float(res)

def calc_rsd(arr, true_val):

  error = 0
  for val in arr:
    error += (val - true_val)**2

  avg_error = np.sqrt(error / len(arr))

  return avg_error/true_val

def calc_rsd_by_series(series, true_val):
    return [calc_rsd(res, true_val) for res in series]

def calc_bias(arr, true_val):

    avg = arr.mean()

    return (avg - true_val)/true_val

def calc_bias_by_series(series, true_val):
    return [calc_bias(res, true_val) for res in series]

def calc_error(arr, true_val):
  avg = arr.mean()
  return round(abs(float(((true_val-avg)/true_val) * 100)), 2)

def error_diagrams(plt, res, res_regul, well):
    # Визуализация
    plt.figure(figsize=(14, 4))

    categories = ['$p_R$', '$r(1)$', '$r(2)$']
    # Настройки отображения
    bar_width = 0.35  # Ширина столбца
    x = np.arange(len(categories))

    # График давления
    ax = plt.subplot(1, 3, 1)

    n = 0
    data1 = [calc_error(res['p_R'][0, 0, n], well.reservoir.p_R),
             calc_error(res['r_1'][0, 0, n], well.reservoir.r_1),
             calc_error(res['r_2'][0, 0, n], well.reservoir.r_2)]
    data2 = [calc_error(res_regul['p_R'][0, 0, n], well.reservoir.p_R),
             calc_error(res_regul['r_1'][0, 0, n], well.reservoir.r_1),
             calc_error(res_regul['r_2'][0, 0, n], well.reservoir.r_2)]
    bars1 = ax.bar(x - bar_width / 2, data1, bar_width, label='Оценка', color='skyblue')
    bars2 = ax.bar(x + bar_width / 2, data2, bar_width, label='Оценка с регул.', color='lightcoral')

    # Добавляем подписи значений над столбцами
    ax.bar_label(bars1, label_type='edge')
    ax.bar_label(bars2, label_type='edge')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    plt.xlabel('Амплитуда шума, 0%')
    plt.ylabel('$\delta,\ \%$')
    plt.grid()
    plt.legend()
    plt.ylim(0, 80)

    # График НЧ шума
    ax = plt.subplot(1, 3, 2)
    n = len(res['p_R'][0, 0]) // 2
    data1 = [calc_error(res['p_R'][0, 0, n], well.reservoir.p_R),
             calc_error(res['r_1'][0, 0, n], well.reservoir.r_1),
             calc_error(res['r_2'][0, 0, n], well.reservoir.r_2)]
    data2 = [calc_error(res_regul['p_R'][0, 0, n], well.reservoir.p_R),
             calc_error(res_regul['r_1'][0, 0, n], well.reservoir.r_1),
             calc_error(res_regul['r_2'][0, 0, n], well.reservoir.r_2)]
    bars1 = ax.bar(x - bar_width / 2, data1, bar_width, label='Оценка', color='skyblue')
    bars2 = ax.bar(x + bar_width / 2, data2, bar_width, label='Оценка с регул.', color='lightcoral')

    # Добавляем подписи значений над столбцами
    ax.bar_label(bars1, label_type='edge')
    ax.bar_label(bars2, label_type='edge')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    plt.xlabel('Амплитуда шума, 5%')
    plt.grid()
    plt.legend()
    plt.ylim(0, 80)

    # График ВЧ шума
    ax = plt.subplot(1, 3, 3)
    n = -1
    data1 = [calc_error(res['p_R'][0, 0, n], well.reservoir.p_R),
             calc_error(res['r_1'][0, 0, n], well.reservoir.r_1),
             calc_error(res['r_2'][0, 0, n], well.reservoir.r_2)]
    data2 = [calc_error(res_regul['p_R'][0, 0, n], well.reservoir.p_R),
             calc_error(res_regul['r_1'][0, 0, n], well.reservoir.r_1),
             calc_error(res_regul['r_2'][0, 0, n], well.reservoir.r_2)]
    bars1 = ax.bar(x - bar_width / 2, data1, bar_width, label='Оценка', color='skyblue')
    bars2 = ax.bar(x + bar_width / 2, data2, bar_width, label='Оценка с регул.', color='lightcoral')

    # Добавляем подписи значений над столбцами
    ax.bar_label(bars1, label_type='edge')
    ax.bar_label(bars2, label_type='edge')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    plt.xlabel('Амплитуда шума, 10%')
    plt.grid()
    plt.legend()
    plt.ylim(0, 80)

    plt.tight_layout()
    return plt
