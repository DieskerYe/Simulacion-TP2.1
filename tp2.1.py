import matplotlib.pyplot as plt
from random import randint
import random as rnd
import math
import random
import secrets
from scipy.stats import chi2
from math import floor
from scipy.stats import ksone
from math import floor
from scipy.stats import ksone
from itertools import tee
from math import sqrt
from scipy.stats import norm

def manejar_excepciones(*args):
    for i, arg in enumerate(args):
        if not (isinstance(arg, int) and arg >= 0):
            raise ValueError("Incorrect argument value on pos {0}.".format(i))

def inicializar_semilla(semilla, max_value):
    if semilla is None:
        semilla = randint(0, max_value)
    return semilla

def metodo_cuadrado_medio(limit, semilla=None, size=10):
    max_value = 10**size-1
    semilla = inicializar_semilla(semilla, max_value)
    manejar_excepciones(limit, semilla, size)
    count = 0
    limite_inferior = int(size / 2)
    limite_superior = limite_inferior + size
    while count < limit:
        semilla = int(str(semilla*semilla).zfill(size*2)[limite_inferior:limite_superior])
        count += 1   # El nro es normalizado al rango [0,1] para simplificar las comparaciones.
        yield semilla / max_value

def random_python(limit, semilla=None):
    if semilla is None:
        random.seed()
    elif semilla >= 0:
        random.seed(semilla)
    manejar_excepciones(limit)
    count = 0
    while count < limit:
        count += 1
        yield random.random()
    
def secrets_python(limit, bits=64):
    manejar_excepciones(limit, bits)
    max_value = 2**bits
    count = 0
    while count < limit:
        count += 1
        yield secrets.randbits(bits) / max_value

def corrida(generador, a=0.05):
    #Generador de corridas de la muestra. Un paso ascendente está representado
    # con un "+", mientras que un paso descendente está representado con un "-".
    corridas = ("+" if a < b else "-" for a, b in pairwise(generador))
    # Contar cantidad de corridas.
    tamaño_muestra = 2 
    cant_corridas = 1
    paso_anterior = next(corridas)
    for paso in corridas:
        tamaño_muestra += 1
        if paso != paso_anterior:
            cant_corridas += 1
        paso_anterior = paso
    # Calcular media y variancia.
    media = (2*tamaño_muestra - 1) / 3
    variancia = (16*tamaño_muestra - 29) / 90
    # Estandarizar distribución normal.
    z0 = (cant_corridas - media) / sqrt(variancia)
    # Obtener valor crítico desde librería ("tabla").
    valor_critico = norm.ppf(q=1-a/2)
    return abs(z0) <= valor_critico

def huecos(generador, base=20, a=0.05):
    # Obtener muestra escalada por un factor base.
    muestra = [min(floor(i*base), base-1) for i in generador]
    # Calcular frecuencia absoluta de cada tamaño de hueco.
    frec_abs = {}
    for nro in range(base):
        posiciones = [index for index, elemento in enumerate(muestra) if elemento == nro]
        for hueco in pairwise(posiciones):
            tamaño_hueco = hueco[1] - hueco[0] - 1
            frec_abs[tamaño_hueco] = frec_abs.get(tamaño_hueco, 0) + 1
    # Calcular frecuencia absoluta acumulada.
    frec_abs_acum = [frec_abs.get(0, 0)]
    for tamaño_hueco in range(1, max(frec_abs.keys())+1):
        frec_abs_acum.append(frec_abs_acum[-1] + frec_abs.get(tamaño_hueco, 0))
    # Calcular frecuencia relativa acumulada.
    frec_rel_acum = []
    for frec in frec_abs_acum:
        frec_rel_acum.append(frec / frec_abs_acum[-1])
    # Calcular diferencia con la frecuecia relativa acumulada esperada.
    n_interv = len(frec_rel_acum)
    frec_esp = [1-(1-1/base)**(x+1) for x in range(n_interv)]
    d = max(abs(frec_rel_acum[i] - frec_esp[i]) for i in range(n_interv))
    # Calcular diferencia de confiabilidad.
    d_conf = ksone.ppf(1-a/2, len(muestra))
    return d < d_conf

def metodo_kolmogorov_smirnov(generador, a=0.05):
    muestra = sorted(generador)
    tamaño_muestra = len(muestra)
    d_max = max(((i+1)/tamaño_muestra - muestra[i]) for i in range(tamaño_muestra))
    d_min = min((muestra[i] - i/tamaño_muestra) for i in range(tamaño_muestra))
    d = max(d_max, abs(d_min))
    valor_critico = ksone.ppf(1-a/2, tamaño_muestra)
    return d < valor_critico

def metodo_chi_cuadrado(generador, n_interv=50, a=0.05):
    muestra = list(generador)
    n_obsv = len(muestra)
    gr_lib = n_interv - 1
    frec_esperada = n_obsv / n_interv
    intervalos = [0] * n_interv
    for elemento in muestra:
        index = floor(elemento * n_interv)
        intervalos[index] += 1
    chi = sum(((frec_observada - frec_esperada)**2 / frec_esperada) \
        for frec_observada in intervalos)
    valor_critico = chi2.isf(q=a, df=gr_lib)
    return chi < valor_critico

def metodo_lcg(limit, semilla=None, m=2**32, a=134775813, c=1):
    max_value = m - 1
    semilla = inicializar_semilla(semilla, max_value)
    manejar_excepciones(limit, semilla, m, a, c)
    count = 0
    while count < limit:
        semilla = (a * semilla + c) % m
        count += 1
        # El nro es normalizado al rango [0,1] para simplificar las comparaciones.
        yield semilla / max_value
def randu(limit, semilla=None):
    return metodo_lcg(limit, semilla=semilla, m=2**31, a=65539, c=0)    

def grafico_dispersion_2d(generador, nros_gen, color='red', **kwargs):
    iterador = generador(nros_gen, **kwargs)
    x = [next(iterador)]
    y = []
    for n in iterador:
        x.append(n)
        y.append(n)
    x.pop()
    ax = plt.axes()
    ax.set_xlabel('Xi')
    ax.set_ylabel('X(i+1)')
    ax.scatter(x, y, color=color)
    plt.show()

def grafico_dispersion_3d(generador, nros_gen, color='red', **kwargs):
    iterador = generador(nros_gen, **kwargs)
    primero = next(iterador)
    segundo = next(iterador)
    x = [primero, segundo]
    y = [segundo]
    z = []
    for n in iterador:
        x.append(n)
        y.append(n)
        z.append(n)
    x.pop()
    x.pop()
    y.pop()
    ax = plt.axes(projection ="3d")
    ax.set_xlabel('Xi')
    ax.set_ylabel('X(i+1)')
    ax.set_zlabel('X(i+2)')
    ax.scatter3D(x, y, z, color=color)
    plt.show()

def generar_reporte(test, gen, n, i, **kwargs):
    positive_results = 0
    for _i in range(i):
        if test(gen(n), **kwargs):
            positive_results += 1
    params = (key + "=" + str(value) for key, value in kwargs.items())
    header = test.__name__ + "(" + ",".join(params) + ")"
    print("{:>40} : {:f}".format(header, positive_results / i))

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

if __name__ == "__main__":
    nros_gen = 1000
    pruebas = 250
    print("~ metodo_lcg ~")
    generar_reporte(metodo_chi_cuadrado, metodo_lcg, nros_gen, pruebas, a=0.02)
    generar_reporte(metodo_kolmogorov_smirnov, metodo_lcg, nros_gen, pruebas, a=0.02)
    generar_reporte(corrida, metodo_lcg, nros_gen, pruebas, a=0.02)
    generar_reporte(huecos, metodo_lcg, nros_gen, pruebas, a=0.02)
    print("~ Middle Square ~")
    generar_reporte(metodo_chi_cuadrado, metodo_cuadrado_medio, nros_gen, pruebas, a=0.02)
    generar_reporte(metodo_kolmogorov_smirnov, metodo_cuadrado_medio, nros_gen, pruebas, a=0.02)
    generar_reporte(corrida, metodo_cuadrado_medio, nros_gen, pruebas, a=0.02)
    generar_reporte(huecos, metodo_cuadrado_medio, nros_gen, pruebas, a=0.02)
    print("~ Random ~")
    generar_reporte(metodo_chi_cuadrado, random_python, nros_gen, pruebas, a=0.02)
    generar_reporte(metodo_kolmogorov_smirnov, random_python, nros_gen, pruebas, a=0.02)
    generar_reporte(corrida, random_python, nros_gen, pruebas, a=0.02)
    generar_reporte(huecos, random_python, nros_gen, pruebas, a=0.02)
    print("~ Secrets ~")
    generar_reporte(metodo_chi_cuadrado, secrets_python, nros_gen, pruebas, a=0.02)
    generar_reporte(metodo_kolmogorov_smirnov, secrets_python, nros_gen, pruebas, a=0.02)
    generar_reporte(corrida, secrets_python, nros_gen, pruebas, a=0.02)
    generar_reporte(huecos, secrets_python, nros_gen, pruebas, a=0.02)

    grafico_dispersion_2d(metodo_lcg, nros_gen)
    grafico_dispersion_2d(metodo_cuadrado_medio, nros_gen)
    grafico_dispersion_2d(random_python, nros_gen)
    grafico_dispersion_2d(secrets_python, nros_gen)

    grafico_dispersion_3d(metodo_lcg, nros_gen)
    grafico_dispersion_3d(metodo_cuadrado_medio, nros_gen)
    grafico_dispersion_3d(random_python, nros_gen)
    grafico_dispersion_3d(secrets_python, nros_gen)