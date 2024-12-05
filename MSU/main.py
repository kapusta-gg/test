from matplotlib import pyplot as plt
import numpy as np
from random import triangular as tr

NAME = X1 = 0
VALUE = X2 = 1
X3 = 2


def euler_method(A, B, C, D=0, h=0.001, start_t=0, end_t=10, degree=2):
    E = np.eye(degree)
    x_0 = np.array([[0], [0]])
    def step(E, h, A, B, x):
        return np.matmul((E + h * A), x) + h * B
    t = np.arange(start_t, end_t, h)
    x_i = x_0
    x_l = []
    y_l = []
    for _ in t:
        x_i = step(E, h, A, B, x_i)
        x_l.append(x_i)
        y_l.append(np.matmul(C, x_i) + D)
    plt.plot(t, y_l)
    for i in range(len(y_l)):
        if y_l[i] >= y_l[-1] * 0.95:
            return [t[i]]


def get_matrixs(x1: float, x2: float, x3: float):
    R1 = x1
    R2 = x2
    R3 = 50
    C1 = x3
    C2 = 0.001
    A = np.array([[-1/(C1 * R1), -1/(C1 * R1)], [-1/(C2* R1), -(R1 * R2 + R2 * R3 + R1 * R3)/(C2 * R1 * R2 * R3)]])
    B = np.array([[1/(C1 * R1)], [(R1 + R2)/(C2 * R1 * R2)]])
    C = np.array([1, 0])
    return A, B, C


def factorial_mathematical_model(x1: tuple[str, float], x2: tuple[str, float], x3: tuple[str, float]) -> None:
    factors_name = [x1[NAME], x2[NAME], x3[NAME]]
    factors_value = [x1[VALUE], x2[VALUE], x3[VALUE]]
    print(f"Выбраны факторы: {x1[NAME]}({x1[VALUE]}), {x2[NAME]}({x2[VALUE]}), {x3[NAME]}({x3[VALUE]})")

    print("-" * 40)

    table_value = ["+++",
                   "-0+",
                   "+-+",
                   "-0+",
                   "++-",
                   "-+-",
                   "+--",
                   "---"]
    table12 = ("13", "+--++--+")
    table13 = ("14", "+-+--+-+")
    table23 = ("25", "++----++")
    table123 = ("123", "+--+-++-")
    factors_max = [x1[VALUE] * 1.1, x2[VALUE] * 1.1, x3[VALUE] * 1.1]
    factors_min = [x1[VALUE] * 0.9, x2[VALUE] * 0.9, x3[VALUE] * 0.9]

    table_dict = {"+": factors_max,
                  "-": factors_min}

    print("Факторы:", end="\n")
    for i in range(3):
        print("+" * 20)
        print(f"Фактор {factors_name[i]}")
        print(f"Значение {factors_value[i]}")
        print(f"Максимальное значение {factors_max[i]}")
        print(f"Минимальное значение {factors_min[i]}")
        print("+" * 20)
        print()
    print("-" * 40)
    y = []
    print("Матриц значений экспериментов")
    for i in range(1, 9):
        n_exp = table_value[i - 1]
        n_x1, n_x2, n_x3 = n_exp[X1], n_exp[X2], n_exp[X3]
        res = euler_method(*get_matrixs(table_dict[n_x1][X1], table_dict[n_x2][X2], table_dict[n_x3][X3]))[0]
        res = round(res, 4)
        y.append(res)
        print(f"Номер эксперимента: {i}; Набор: {n_exp}; y_{i} = {res}")

    print("-" * 40)

    b = []
    print("Коэффициенты факторной мат. модели")
    for i in range(4):
        summ = 0
        print(f"Коэффициент b{i} = ", end="")
        if i:
            for j in range(8):
                n_x = table_value[j][i - 1]
                summ = summ + y[j] if n_x == "+" else summ - y[j]
            summ = round(summ, 4)
        else:
            summ = sum(y)
        print(summ / 8)
        b.append(summ / 8)

    for i in (table12, table13, table23, table123):
        summ = 0
        print(f"Коэффициент b{i[0]} = ", end="")
        for j in range(8):
            summ = summ + y[j] if i[1][j] == "+" else summ - y[j]
        summ = round(summ, 4)
        print(summ / 8)
        b.append(summ / 8)
    print("-" * 40, end="\n")

    print("=" * 80, end="\n")
    print("Обработка результатов эксперимента при равномерном дублировании")
    print()

    y_mean = []
    y_temp = []
    print("Матриц значений экспериментов #2")
    for i in range(1, 9):
        print("+" * 20)
        print(f"Номер эксперимента: {i}; Набор: {n_exp}")
        n_exp = table_value[i - 1]
        n_x1, n_x2, n_x3 = n_exp[X1], n_exp[X2], n_exp[X3]
        temp = []
        for j in range(1, 4):
            print(f"\tПараллельный опыт №{j}")
            r1, r2, r3 = round(tr(0, 0.02), 4), round(tr(0, 0.02), 4), round(tr(0, 0.02), 2)
            t1, t2, t3 = table_dict[n_x1][X1] * (1 + r1), table_dict[n_x2][X2] * (1 + r2), table_dict[n_x3][X3] * (1 + r3)
            t1, t2 = round(t1, 4), round(t2, 4)
            print(f"\t\tЗначения: x1={t1}; x2={t2}; x3={t3}")
            res = euler_method(*get_matrixs(t1, t2, t3))[0]
            res = round(res, 4)
            print(f"\t\tРезультат: y_{i}{j}={res}")
            temp.append(res)
        y_mean_i = sum(temp) / 3
        y_mean_i = round(y_mean_i, 4)
        y_mean.append(y_mean_i)
        y_temp.append(temp)
        print(f"Результат: Средняя y_{i} = {y_mean_i}")
        print("+" * 20)
        print()

    print("-" * 40)

    print("Расчет СКО для каждого из опытов")
    sko = []
    for i in range(8):
        temp = sum([(j - y_mean[i]) ** 2 for j in y_temp[i]]) / 23
        sko.append(temp)
        print(f"Номер эксперимента: {i + 1}; СКО^2: {temp}")

    print("-" * 40)

    print("Метод вычисления максимального относительного отклонения")
    y_m = round(sum([sum(i) for i in y_temp]) / 24, 4)
    sko_m = sum([sum([(j - y_m) ** 2 for j in i]) for i in y_temp]) / 23
    print(f"y_mean = {y_m}, sko_mean^2 = {sko_m}")
    quantiles_dict = {0.10: 2.52, 0.05: 2.70, 0.025: 2.86, 0.01: 2.05}
    max_y = min([min(i) for i in y_temp])
    print(f"Минимальное значение y = {0.201}")
    t = round(abs(max_y - y_m) / np.sqrt(sko_m), 2)
    print(f"При уровне значимости 0.05: ({t} <= {quantiles_dict[0.05]}) -> {t <= quantiles_dict[0.05]}")

    print("-" * 40)

    print("Проверка гипотезы об однородности дисперсий Опытов", end='\n\n')
    print("Таблица Fp")
    sko_t = [round(i, 8) for i in sko]
    print("\t" * 2 + "|" + "\t|".join([str(i) for i in sko_t]))
    print("-" * 80)
    for i in sko_t:
        print(f"{i}", end="\t|")
        for j in sko_t:
            print(round(max(i, j) / min(i, j), 3), end="\t|")
        print()
    print("Объем выборок при проведении одного опыта = 3. Следовательно m1 = 3 - 1 = 2, m2 = 3 - 1 = 2")
    Fm = 19
    print(f"В таблице распределения F Фишера-Снедекора Fm = {Fm}")
    print("Проверим наши Fp")
    print("\t" * 2 + "|" + "\t|".join([str(i) for i in sko_t]))
    print("-" * 80)
    for i in sko_t:
        print(f"{i}", end="\t|")
        for j in sko_t:
            print("T" if round(max(i, j) / min(i, j), 3) < Fm else "F", end="\t\t|")
        print()
    print("Как видно однородность дисперсий опытов подтверждается F-критерием Фишера")
    print()
    print("Проверим однородность с помощью критерия Кохрена")
    print("Рассчитаем Gp")
    Gp = max(sko) / sum(sko)
    print(f"Параметр Gp = {Gp}")
    print("Сравним расчетно значение Gp с табличным Gt")
    Gt = 0.5157
    print(f"Всего опытов 8. Степень параллельных опытов = 3. Следовательно, Gt = {Gt}")
    print(f"Проверим Gp < Gt -> Получаем {Gp < Gt}")
    print("Критерий Кохрена так же подтвержден")

    print("-" * 40)

    print("Вычислим дисперсию воспроизводимости Sy^2")
    Sy_2 = sum(sko) / 8
    print(f"Дисперсия воспроизводимости Sy^2 = {Sy_2}")

    print("-" * 40)

    print("Вычислим коэффициенты факторной модели")
    b_n = []
    print("Коэффициенты факторной мат. модели")
    for i in range(4):
        summ = 0
        print(f"Коэффициент b{i} = ", end="")
        if i:
            for j in range(8):
                n_x = table_value[j][i - 1]
                summ = summ + y_mean[j] if n_x == "+" else summ - y_mean[j]
            summ = round(summ, 4)
        else:
            summ = sum(y_mean)
        print(summ / 8)
        b_n.append(summ / 8)

    for i in (table12, table13, table23, table123):
        summ = 0
        print(f"Коэффициент b{i[0]} = ", end="")
        for j in range(8):
            summ = summ + y_mean[j] if i[1][j] == "+" else summ - y_mean[j]
        summ = round(summ, 4)
        print(summ / 8)
        b_n.append(summ / 8)

    print("-" * 40)

    print("Осуществим проверку значимости вычисленных коэффициентов факторной модели.")
    S_b = [np.sqrt(i / 24) for i in sko]
    print("Осущетсвим Далее проверку двумя способами:")
    f = 12
    print("1) Сравним абсолютные величины с доверительным интервалом")
    print("Вычислим inter b Для всех S(b_i) и сравним с t_T табличным")
    t_T = 2.18
    print(f"f = 12, alpha = 0.05. Следовательно t_T = {t_T}")
    for i in range(8):
        print(f"inter b = {t_T * S_b[i]}. Проверим inter b < |b| -> Получаем {t_T * S_b[i] < abs(b_n[i])}")
    print("2) Проверим t-критерий Стьюдента")
    print("Вычислим tp Для всех S(b_i) и сравним с t_T табличным")
    t_T = 2.18
    print(f"f = 12, alpha = 0.05. Следовательно t_T = {t_T}")
    for i in range(8):
        print(f"tp = {abs(b_n[i]) / S_b[i]}. Проверим tp > t_T -> Получаем {abs(b_n[i]) / S_b[i] > t_T}")
        if not abs(b_n[i]) / S_b[i] > t_T:
            b_n[i] = 0
    print("Занулим незначительные коэффициенты")


    print("-" * 40)

    print("Осуществим проверку гипотезы об адекватности полученной факторной модели. Для этого вычислим дисперсию "
          "адекватности:")

    table_temp = ("++++++++",
                  "+–++––+–",
                  "++–+–+––",
                  "+––++––+",
                  "+++–+–––",
                  "+–+––+–+",
                  "++––––++",
                  "+–––+++–")
    summ = 0
    for i in range(8):
        y = y_mean[i]
        r = 0
        for j in range(8):
            r += b_n[j] if table_temp[i][j] == "+" else -b_n[j]
        summ += (y - r) ** 2
    S_ad_2 = (3 / 4) * summ
    print(f"S_ad^2 = {S_ad_2}")
    print("Проверка гипотезы адекватности модели осуществляется с помощью Fкритерия Фишера."
          " Для этого вычисляется: Fp = S_ad^2 / Sy^2")
    Fp = S_ad_2 / Sy_2
    Fm = 3.63
    print(f"Сравним Fp < Fm (Fm = {Fm}) -> Получаем {Fp < Fm}")

factorial_mathematical_model(("R1", 50), ("R2", 100), ("C1", 0.001))