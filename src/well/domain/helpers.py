"""Модуль вспомогательных функций"""


def transfer_function(
    current_val: float, target: float, tf_time: float = 1, modeling_dt: float = 1
):
    """Апериодическое звено первого порякда"""
    return current_val + ((modeling_dt / tf_time) * (target - current_val))


FILTER_MAP = {"aperiodic": transfer_function}


def filter_data(
    data: list, filter_type: str, modeling_dt: float, tf_time: float
) -> list:
    """Фильтрация данных"""

    filter_func = FILTER_MAP.get(filter_type)

    if not filter_func:
        return data

    res = []

    current_val = 0

    for k, val in enumerate(data):

        if k == 0:
            current_val = val
        else:
            current_val = filter_func(current_val, val, tf_time, modeling_dt)

        res.append(current_val)

    return res
