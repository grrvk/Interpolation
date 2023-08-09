# This is a sample Python script.
import functions as f

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
polinomial = f.np.zeros((f.N, f.N))
multipliers = f.identity_row_creation(f.N)
diffs = f.calc_divs()

polinomial = f.fill_polinomial_matrix(polinomial, diffs)

polinomial_array = f.create_polinomial_array(polinomial)

f.print_polinomial(polinomial_array)
print("\n\n---------------------------------------------------------------------\n")
'''---------------------------------------------------------'''

polinomial_ch = f.np.zeros((f.N, f.N))
diffs_ch = f.calc_divs_ch()

polinomial_ch = f.fill_polinomial_ch_matrix(polinomial_ch, diffs_ch)

polinomial_array_ch = f.create_polinomial_array(polinomial_ch)
f.print_polinomial(polinomial_array_ch)
print("\n\n---------------------------------------------------------------------\n")
'''---------------------------------------------------------'''
x, y = f.f_results()
x_E, y_E = f.E_results(polinomial_array)
x_T, y_T = f.T_results(polinomial_array_ch)
f.full_plot(x, y, x_E, y_E, x_T, y_T)
f.decs(x, y, y_E, y_T)
'''---------------------------------------------------------'''

f_arr = f.f_array()
A = f.A_matrix()
H = f.H_matrix()
b = f.mult_hf(H, f_arr)

Augmented_matrix = f.np.zeros((f.N-2, f.N-1))
for i in range(f.N-2):
    for j in range(f.N-2):
        Augmented_matrix[i][j] = A[i][j]
    Augmented_matrix[i][f.N-2] = b[i][0]

P = f.P_identity()
M = f.np.zeros((f.N-2, f.N-2))

for i in range(f.N-2):
    P = f.P_identity()
    M = f.P_identity()
    row_index = 0

    row_index = f.findMax(Augmented_matrix, i, row_index, 0)
    P[[i, row_index]] = P[[row_index, i]]

    Augmented_matrix = f.mult(Augmented_matrix, P)

    M = f.change_M_matrix(M, Augmented_matrix, i)

    Augmented_matrix = f.mult(Augmented_matrix, M)

'''for i in range(f.N-2):
    row = ""
    for j in range(f.N-2):
        row += str(Augmented_matrix[i][j]) + " "
    row += str(Augmented_matrix[i][f.N-2]) + " "
    print(row)'''

X = f.result_fill(Augmented_matrix)
m_full = f.np.zeros((1, f.N))
m_full[0][0] = m_full[0][f.N-1] = 0

for i in range(1, f.N-1):
    m_full[0][i] = X[i-1][0]

g_x = f.g_x(m_full)
spl = f.print_g_matrix(g_x, x_E)
print("\n")

for i in range(f.N-1):
    f.print_spl_row(spl[i], i)

y_spl = f.np.zeros((1, pow(f.N, 3)))

for i in range(pow(f.N, 3)):
    y_spl[0][i] = f.get_spl_y(spl, x[0][i])


f.full_spl_plot(x, y, y_spl)
f.spl_dec(x, y, y_spl)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
