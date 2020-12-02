
from sympy import symbols, diff, integrate


class matriz():

    def __init__(self, grauElem, nelem, L, length):
        import numpy as np
        x = symbols('x')

        self.length = length
        self.npoints = grauElem * nelem + 1

        self.N_x = []  # lista com as funcoes N(x) pra cada elemento
        self.deriv = []  # lista com as derivadas de N(x) pra cada elemento
        if grauElem == 1:
            # Conjunto de coord y dos pontos da N(x) pra cada elemento
            Y = [[1, 0], [0, 1]]
            for elem in range(grauElem + 1):
                for N in range(grauElem + 1):
                    y = Y[N]
                    # Conjunto de coord x dos pontos da N(x)
                    x_p = [elem/nelem, elem/nelem + self.length]
                    # coeficientes da funcao N(x)
                    coefs_N = np.polyfit(x_p, y, grauElem)
                    #  Preenche a lista de funções N(x)
                    self.N_x.append(coefs_N[0] * x + coefs_N[1])
                    #  Preenche a lista de derivadas de N(x)
                    self.deriv.append(diff(self.N_x[N]))

        elif grauElem == 2:
            # Conjunto de coord y dos pontos da N(x) pra cada elemento
            Y = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            for elem in range(grauElem + 1):
                for N in range(grauElem + 1):
                    y = Y[N]
                    # Conjunto de coord x dos pontos da N(x)
                    x_p = [elem/nelem, elem/nelem +
                           self.length/2, elem/nelem + self.length]
                    # coeficientes da funcao N(x)
                    coefs_N = np.polyfit(x_p, y, grauElem)
                    #  Preenche a lista de funções N(x)
                    self.N_x.append(coefs_N[0] * x **
                                    2 + coefs_N[1] * x + coefs_N[2])
                    #  Preenche a lista de derivadas de N(x)
                    self.deriv.append(diff(self.N_x[N]))

        elif grauElem == 3:
            # Conjunto de coord y dos pontos da N(x) pra cada elemento
            Y = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
            for elem in range(grauElem + 1):
                for N in range(grauElem + 1):
                    y = Y[N]
                    # Conjunto de coord x dos pontos da N(x)
                    x_p = [elem/nelem, elem/nelem + self.length/3, 
                           elem/nelem + 2*self.length/3,
                           elem/nelem + self.length]  
                    # coeficientes da funcao N(x)
                    coefs_N = np.polyfit(x_p, y, grauElem)
                    #  Preenche a lista de funções N(x)
                    self.N_x.append(
                        coefs_N[0] * x ** 3 + coefs_N[1] * x ** 2 
                        + coefs_N[2] * x + coefs_N[3])
                    #  Preenche a lista de derivadas de N(x)
                    self.deriv.append(diff(self.N_x[N]))

        # print('N(x)_list =', self.N_x)
        # print('derivatives =', self.deriv)
        self.grauElem = grauElem

    def matrizk(self):
        import numpy as np
        x = symbols('x')
        # Montando a matriz K (são iguais para todos os elementos)
        self.k = np.zeros(
            (self.grauElem + 1, self.grauElem + 1), dtype='float')
        limit_f = self.length
        limit_i = 0
        for i in range(self.grauElem + 1):
            for j in range(self.grauElem + 1):
                f = integrate(self.deriv[i] * self.deriv[j], x)
                # Aplico os valores de x na integral definida
                result = f.subs(x, limit_f) - f.subs(x, limit_i)
                self.k[i, j] = result
        return self.k

    def matrizm(self):
        import numpy as np
        x = symbols('x')
        self.m = np.zeros(
            (self.grauElem + 1, self.grauElem + 1), dtype='float')
        limit_f = self.length
        limit_i = 0
        for i in range(self.grauElem + 1):
            for j in range(self.grauElem + 1):
                f = integrate(self.N_x[i] * self.N_x[j], x)
                # Aplico os valores de x na integral definida
                result = f.subs(x, limit_f) - f.subs(x, limit_i)
                self.m[i, j] = result
        return self.m

    def matrizf(self, fonte):
        import numpy as np
        x = symbols('x')
        self.f = np.zeros((self.grauElem + 1), dtype='float')
        limit_f = self.length
        limit_i = 0
        for j in range(self.grauElem + 1):
            f = integrate(fonte * self.N_x[j], x)
            # Aplico os valores de x na integral definida
            result = f.subs(x, limit_f) - f.subs(x, limit_i)
            self.f[j] = result
        return self.f

    def matrizg(self):
        import numpy as np
        x = symbols('x')
        self.g = np.zeros(
            (self.grauElem + 1, self.grauElem + 1), dtype='float')
        limit_f = self.length
        limit_i = 0
        for i in range(self.grauElem + 1):
            for j in range(self.grauElem + 1):
                f = integrate(self.N_x[j] * self.deriv[i], x)
                # Aplico os valores de x na integral definida
                result = f.subs(x, limit_f) - f.subs(x, limit_i)
                self.g[i, j] = result
        return self.g
