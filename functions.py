import math
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

np.set_printoptions(threshold=sys.maxsize)

N = int(10)  # 11 - має бути кількість вузлів
A = float(-4)  # -4
B = float(5)  # 5
h = float((B - A) / (N - 1))  # 0.66


def count_x(index):
    if index == 0:
        return A
    return A + index * h


def factorial(i):
    if i == 1:
        return i
    return i * factorial(i - 1)


def count_y(x):
    if A <= x <= -2:
        return x + 2
    elif -2 <= x <= 0:
        return -pow((x + 1), 2) + 1
    elif 0 <= x <= 2:
        return pow(x, 3)
    else:
        return 12 - 2 * x


def div_diff(power, index):
    if power == 1:
        x_i1 = float(count_x(index + 1))
        x_i = float(count_x(index))
        return count_y(x_i1) - count_y(x_i)  # count_y
    else:
        return div_diff(power - 1, index + 1) - div_diff(power - 1, index)


def diff_matrix():
    creation = np.zeros((N + 1, N + 1))
    for i in range(N + 1):
        creation[i][0] = count_y(count_x(i))
        for j in range(1, N + 1):
            creation[i][j] = div_diff(j, i)


def identity_row_creation(n):
    creation = np.ones((1, n))
    return creation


def print_matrix(result):
    for i in range(N):
        row = ""
        for j in range(N - i):
            row += str(result[i][j]) + " "
        print(row)


def multiplyTwoPolynomials(A, B, m, n):
    productPolynomial = np.zeros((1, m + n - 1))
    for i in range(m + n - 1):
        productPolynomial[0][i] = 0
    for i in range(m):
        for j in range(n):
            productPolynomial[0][i + j] += A[0][i] * B[0][j]

    return productPolynomial


def copy_row(A, B, row, mult):
    for j in range(row + 1):
        B[row][j] = mult * A[0][row - j]
    return B


def create_polinomial_array(A):
    creation = np.zeros((1, N))

    for j in range(N):
        summ = 0
        for i in range(N):
            if i % 2 != 0:
                summ -= A[i][j]
            else:
                summ += A[i][j]
        creation[0][j] = summ
    return creation


def print_polinomial(polinomial_array):
    print("P(x) = " + str(polinomial_array[0][0]), end="")
    for i in range(1, N):
        if i == int(N / 2) + 1:
            print("\n", end="")
        if i % 2 != 0:
            if polinomial_array[0][i] < 0:
                print(" + " + str(-polinomial_array[0][i]) + "*x^" + str(i), end="")
            if polinomial_array[0][i] > 0:
                print(" - " + str(polinomial_array[0][i]) + "*x^" + str(i), end="")
        else:
            if polinomial_array[0][i] < 0:
                print(" - " + str(-polinomial_array[0][i]) + "*x^" + str(i), end="")
            if polinomial_array[0][i] > 0:
                print(" + " + str(polinomial_array[0][i]) + "*x^" + str(i), end="")


def calc_divs():
    creation = np.zeros((1, N))
    sep_div = np.zeros((N, N))

    for j in range(N):
        sep_div[j][0] = count_y(count_x(j))

    for j in range(1, N):
        for i in range(N - j):
            sep_div[i][j] = ((sep_div[i + 1][j - 1] - sep_div[i][j - 1]) / (count_x(i + j) - count_x(i)));

    for i in range(N):
        creation[0][i] = sep_div[0][i]

    return creation


def fill_polinomial_matrix(polinomial, multipliers):
    equation = "P(x) = " + str(count_y(A))
    polinomial[0][0] = count_y(A)

    for power in range(1, N):
        if (power == N - 6 or power == N - 4 or power == N - 3 or power == N - 2 or power == N - 1):
            equation += "\n"
        equation += " + " + str(multipliers[0][power])
        in_brackets = identity_row_creation(power)

        for j in range(power):
            in_brackets[0][j] = count_x(j)
            if count_x(j) + j > 0:
                equation += "(x - " + str(count_x(j)) + ")"
            elif count_x(j) + j < 0:
                equation += "(x + " + str(-(count_x(j))) + ")"
            else:
                equation += "x"

        if power == 1:
            polinomial[power][0] = in_brackets[0] * multipliers[0][1]
            polinomial[power][1] = multipliers[0][1]
        else:
            max_len = 2
            first = identity_row_creation(power + 1)
            first[0][1] = in_brackets[0][0]
            for i in range(power - 1):
                second = identity_row_creation(2)
                second[0][1] = in_brackets[0][i + 1]
                first = multiplyTwoPolynomials(first, second, max_len, 2)
                max_len = max_len + 1
            copy_row(first, polinomial, power, multipliers[0][power])
    equation += "\n"
    print(equation)
    return polinomial


def count_E_polinomial(polinomial_array, x_i):
    summ = float(0)
    for i in range(N):
        if (i % 2 != 0):
            if polinomial_array[0][i] * pow(x_i, i) < 0:
                summ += -polinomial_array[0][i] * pow(x_i, i)
            if polinomial_array[0][i] * pow(x_i, i) > 0:
                summ -= polinomial_array[0][i] * pow(x_i, i)
        else:
            if polinomial_array[0][i] * pow(x_i, i) < 0:
                summ -= -polinomial_array[0][i] * pow(x_i, i)
            if polinomial_array[0][i] * pow(x_i, i) > 0:
                summ += polinomial_array[0][i] * pow(x_i, i)
    return summ


''' --------------------------------------------------------------- '''


def ch_x(k):
    return (A + B) / 2 + ((B - A) / 2) * math.cos(((2 * k + 1) * math.pi) / (2 * N))


def div_diff_ch(power, index):
    if (power == 1):
        x_i1 = ch_x(index + 1)
        x_i = ch_x(index)
        return count_y(x_i1) - count_y(x_i)
    else:
        return div_diff_ch(power - 1, index + 1) - div_diff_ch(power - 1, index)


def calc_divs_ch():
    creation = np.zeros((1, N))
    sep_div = np.zeros((N, N))
    for j in range(N):
        sep_div[j][0] = count_y(ch_x(j))

    for j in range(1, N):
        for i in range(N - j):
            sep_div[i][j] = ((sep_div[i + 1][j - 1] - sep_div[i][j - 1]) / (ch_x(i + j) - ch_x(i)));

    for i in range(N):
        creation[0][i] = sep_div[0][i]

    return creation


def fill_polinomial_ch_matrix(polinomial, multipliers):
    equation = "P(x) = " + str(count_y(ch_x(0)))
    polinomial[0][0] = count_y(ch_x(0))

    for power in range(1, N):
        if (
                power == N - 7 or power == N - 6 or power == N - 5 or power == N - 4 or power == N - 3 or power == N - 2 or power == N - 1):
            equation += "\n"
        equation += " + " + str(multipliers[0][power])

        in_brackets = identity_row_creation(power)

        for j in range(power):
            in_brackets[0][j] = ch_x(j)
            if ch_x(j) > 0:
                equation += "(x - " + str(ch_x(j)) + ")"
            elif ch_x(j) < 0:
                equation += "(x + " + str(-(ch_x(j))) + ")"
            else:
                equation += "x"

        if power == 1:
            polinomial[power][0] = in_brackets[0] * multipliers[0][1]
            polinomial[power][1] = multipliers[0][1]
        else:
            max_len = 2
            first = identity_row_creation(power + 1)
            first[0][1] = in_brackets[0][0]
            for i in range(power - 1):
                second = identity_row_creation(2)
                second[0][1] = in_brackets[0][i + 1]

                first = multiplyTwoPolynomials(first, second, max_len, 2)
                max_len = max_len + 1
            copy_row(first, polinomial, power, multipliers[0][power])
    equation += "\n"
    print(equation)
    return polinomial


def count_ch_polinomial(polinomial_array_ch, x_i):
    summ = float(0)
    for i in range(N):
        if (i % 2 != 0):
            if polinomial_array_ch[0][i] * pow(x_i, i) < 0:
                summ += -polinomial_array_ch[0][i] * pow(x_i, i)
            if polinomial_array_ch[0][i] * pow(x_i, i) > 0:
                summ -= polinomial_array_ch[0][i] * pow(x_i, i)
        else:
            if polinomial_array_ch[0][i] * pow(x_i, i) < 0:
                summ -= -polinomial_array_ch[0][i] * pow(x_i, i)
            if polinomial_array_ch[0][i] * pow(x_i, i) > 0:
                summ += polinomial_array_ch[0][i] * pow(x_i, i)
    return summ


''' --------------------------------------------------------------- '''


def f_results():
    x = identity_row_creation(pow(N, 3))
    y = identity_row_creation(pow(N, 3))
    for i in range(pow(N, 3)):
        x[0][i] = A + ((B - A) / (pow(N, 3) - 1)) * i
    for i in range(pow(N, 3)):
        y[0][i] = count_y(x[0][i])

    return x, y


def E_results(polinomial_array):
    x = identity_row_creation(pow(N, 3))
    y = identity_row_creation(pow(N, 3))
    for i in range(pow(N, 3)):
        x[0][i] = A + ((B - A) / (pow(N, 3) - 1)) * i
    for i in range(pow(N, 3)):
        y[0][i] = count_E_polinomial(polinomial_array, x[0][i])

    return x, y


def T_results(polinomial_array):
    x = identity_row_creation(pow(N, 3))
    y = identity_row_creation(pow(N, 3))
    for i in range(pow(N, 3)):
        x[0][i] = A + ((B - A) / (pow(N, 3) - 1)) * i
    for i in range(pow(N, 3)):
        y[0][i] = count_ch_polinomial(polinomial_array, x[0][i])

    return x, y


def full_plot(x, y, x_E, y_E, x_T, y_T):
    plt.suptitle('Усі графіки')
    plt.plot(x[0], y[0], label="f(x)")
    plt.plot(x_E[0], y_E[0], label="Рівновіддалені вузли")
    plt.plot(x_T[0], y_T[0], label="Чебишовські вузли")
    plt.legend()
    plt.show()


def decs(x, y, y_E, y_T):
    plt.suptitle("Різниці")
    plt.plot(x[0], y[0] - y_E[0], label="f(x) - E", color = 'orange')
    plt.plot(x[0], y[0] - y_T[0], label="f(x) - T", color = 'green')
    plt.legend()
    plt.show()



'''--------------------------------------'''


def mult_hf(H, f_x):
    res = np.zeros((N - 2, 1))

    for i in range(N-2):
        for j in range(N):
            res[i][0] += H[i][j] * f_x[j][0]

    return res

def P_identity():
    P = np.zeros((N-2, N-2))
    for i in range(N-2):
        P[i][i] = 1
    return P

def findMax(A, col_index, row_index, max):
    for i in range(col_index, N-2):
        if (abs(A[i][col_index]) > max):
            max = A[i][col_index]
            row_index = i
    return row_index


def mult(Augmented_matrix, P):
    res = np.zeros((N - 2, N-1))

    for i in range(N-2):
        for j in range(N-1):
            for k in range(N-2):
                res[i][j] += P[i][k] * Augmented_matrix[k][j]

    return res


def change_M_matrix(M, A, col_num):
    for i in range (col_num, N-2):
        num = A[i][col_num]
        mid_num = A[col_num][col_num]
        if (i != col_num):
            if (num==0):
                M[i][col_num] = 0
            else:
                if (mid_num!=0):
                    M[i][col_num] = -num/mid_num
        else:
            if (mid_num!=0):
                M[i][col_num] = 1/mid_num
    return M


def result_fill(Augmented_matrix):
    X = np.zeros((N-2, 1))
    for i in range(N-3, -1, -1):
        main = Augmented_matrix[i][N-2]
        X[i] = main
        for m in range(i+1, N-2):
            if (Augmented_matrix[i][m]!=0):
                X[i] =  main - Augmented_matrix[i][m]*X[m]
                main -= Augmented_matrix[i][m]*X[m]

    return X



'''--------------------------------------'''


def f_array():
    f = np.zeros((N, 1))
    for j in range(N):
        f[j][0] = count_y(count_x(j))

    return f


def A_matrix():
    A = np.zeros((N - 2, N - 2))
    for i in range(1, N - 1):
        A[i - 1][i - 1] = float(2 * h / 3)
        if i - 2 >= 0:
            A[i - 1][i - 2] = float(h / 6)
        if i < N - 2:
            A[i - 1][i] = float(h / 6)

    return A


def H_matrix():
    H = np.zeros((N - 2, N))

    for i in range(1, N - 1):
        H[i - 1][i - 1] = float(1 / h)
        H[i - 1][i] = float(-2 / h)
        H[i - 1][i + 1] = float(1 / h)
    return H


def g_x(m_full):
    g_matrix = np.zeros((N - 1, 4))

    for i in range(1, N):
        g_matrix[i - 1][0] = float(m_full[0][i - 1] / (6 * h))
        g_matrix[i - 1][1] = float(m_full[0][i] / (6 * h))
        g_matrix[i - 1][2] = float((6 * count_y(count_x(i - 1)) - m_full[0][i - 1] * h * h) / (6 * h))
        g_matrix[i - 1][3] = float((6 * count_y(count_x(i)) - m_full[0][i] * h * h) / (6 * h))

    return g_matrix


def print_g_matrix(g_matrix, x_E):
    spl_matrix = np.zeros((N - 1, 4))

    for i in range(1, N):
        row = "g(x) = " + str(g_matrix[i - 1][0]) + " *(" + str(count_x(i)) + " - x)^3"
        if g_matrix[i - 1][1] < 0:
            if x_E[0][i - 1] > 0:
                row += " - " + str(-g_matrix[i - 1][1]) + " *(x - " + str(count_x(i - 1)) + ")^3"
            elif x_E[0][i - 1] < 0:
                row += " - " + str(-g_matrix[i - 1][1]) + " *(x + " + str(-count_x(i - 1)) + ")^3"
            else:
                row += " - " + str(-g_matrix[i - 1][1]) + " *x^3"
        if g_matrix[i - 1][1] > 0:
            if count_x(i - 1) > 0:
                row += " + " + str(g_matrix[i - 1][1]) + " *(x - " + str(count_x(i - 1)) + ")^3"
            elif count_x(i - 1) < 0:
                row += " + " + str(g_matrix[i - 1][1]) + " *(x + " + str(-count_x(i - 1)) + ")^3"
            else:
                row += " + " + str(g_matrix[i - 1][1]) + " *x^3"
        if g_matrix[i - 1][2] < 0:
            row += " - " + str(-g_matrix[i - 1][2]) + " *(" + str(count_x(i)) + " - x)"
        if g_matrix[i - 1][2] > 0:
            row += " + " + str(g_matrix[i - 1][2]) + " *(" + str(count_x(i)) + " - x)"
        if g_matrix[i - 1][3] < 0:
            if count_x(i - 1) > 0:
                row += " - " + str(-g_matrix[i - 1][3]) + " *(x - " + str(count_x(i - 1)) + ")"
            elif count_x(i - 1) < 0:
                row += " - " + str(-g_matrix[i - 1][3]) + " *(x + " + str(-count_x(i - 1)) + ")"
            else:
                row += " - " + str(-g_matrix[i - 1][3]) + " *x"
        if g_matrix[i - 1][3] > 0:
            if count_x(i - 1) > 0:
                row += " + " + str(g_matrix[i - 1][3]) + " *(x - " + str(count_x(i - 1)) + ")"
            elif count_x(i - 1) < 0:
                row += " + " + str(g_matrix[i - 1][3]) + " *(x + " + str(-count_x(i - 1)) + ")"
            else:
                row += " + " + str(g_matrix[i - 1][3]) + " *x"

        spl_matrix[i - 1][0] = g_matrix[i - 1][1] - g_matrix[i - 1][0]
        spl_matrix[i - 1][1] = g_matrix[i - 1][0] * 3 * count_x(i) - g_matrix[i - 1][1] * 3 * count_x(i - 1)
        spl_matrix[i - 1][2] = 3 * g_matrix[i - 1][1] * count_x(i - 1) * count_x(i - 1) - g_matrix[i - 1][
            0] * 3 * count_x(
            i) * count_x(i) - g_matrix[i - 1][2] + g_matrix[i - 1][3]
        spl_matrix[i - 1][3] = g_matrix[i - 1][0] * count_x(i) * count_x(i) * count_x(i) - g_matrix[i - 1][1] * count_x(
            i - 1) * count_x(i - 1) * count_x(i - 1) + g_matrix[i - 1][2] * count_x(i) - g_matrix[i - 1][3] * count_x(
            i - 1)

        row += " , x є [" + str(count_x(i - 1)) + ", " + str(count_x(i)) + "]"
        print(row)
    return spl_matrix


def print_spl_row(spl_matrix, j):
    row = "g(x) = " + str(spl_matrix[0]) + "*x^3"
    for i in range(1, 4):
        if spl_matrix[i] < 0:
            row += " - " + str(-spl_matrix[i]) + "*x^" + str(3 - i)
        if spl_matrix[i] > 0:
            row += " + " + str(spl_matrix[i]) + "*x^" + str(3 - i)
    row += " , x є [" + str(count_x(j)) + ", " + str(count_x(j + 1)) + "]"
    print(row)


def get_spl_y(spl_matrix, x):
    i = 0
    while x >= count_x(i):
        i = i + 1
    if (x == B):
        i = i - 1

    return (spl_matrix[i - 1][0] * pow(x, 3)) + (spl_matrix[i - 1][1] * pow(x, 2)) + (spl_matrix[i - 1][2] * x) + \
        spl_matrix[i - 1][3]


def full_spl_plot(x, y, y_E):
    plt.suptitle('Усі графіки')
    plt.plot(x[0], y[0], label="f(x)")
    plt.plot(x[0], y_E[0], label="s(x)")
    plt.legend()
    plt.show()


def spl_dec(x, y, y_T):
    plt.suptitle("f(x) - s(x)")
    plt.plot(x[0], y[0] - y_T[0])
    plt.show()
