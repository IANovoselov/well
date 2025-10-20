"""Модуль вспомогательных функций"""
from idlelib.editor import darwin


def transfer_function(current_val: float, target: float, tf_time: float = 1, dt: float = 1):
    return current_val + ((dt / tf_time) * (target - current_val))

FILTER_MAP = {'aperiodic': transfer_function}


def filter_data(data: list, filter_type: str, dt: float, tf_time: float) -> list:
    """Фильтрация данных"""

    filter_func = FILTER_MAP.get(filter_type)

    if not filter_func:
        return data

    res = []

    for k, val in enumerate(data):

        if k == 0:
            current_val = val
        else:
            current_val = filter_func(current_val, val, tf_time, dt)

        res.append(current_val)

    return res


