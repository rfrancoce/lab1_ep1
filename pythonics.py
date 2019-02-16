import numpy as np
#Iniciar la semilla para garantizar que los datos serán iguales cada que se corra el algoritmo
np.random.seed(204)
todos_turnos = []
#Definir cuántas veces se corre la simulación
muestras = 600
for x in range(muestras) :
    #Comenzar el turno sin monedas
    monedas = 0
    turno_aleatorio = [0]
    for x in range(100) :
        dado = np.random.randint(1,7)
        if dado <= 3 :
            monedas = max(0,monedas - 1)
        elif dado < 6 :
            monedas = monedas + 1
        else :
            monedas = monedas + np.random.randint(1,7)
        #Registrar cuantas monedas tengo al final de cada tirada
        turno_aleatorio.append(monedas)
    #Guardar los resultados del turno
    todos_turnos.append(turno_aleatorio)

#Formatear el arreglo como numpy array 
np_todos_turnos = np.array(todos_turnos) 

#Trasponer filas por columnas para adaptar a la gráfica 
np_todos_turnos_t = np.transpose(np_todos_turnos) 

#Sacar la última fila - resultado final de todos los turnos 
ultimos = np_todos_turnos_t[-1,:] 

#Calcular probabilidad de ganar contando los valores del vector
#mayores o iguales a 50 y dividiendo por el número de turnos 
print('La probabilidad de ganar el juego es de ' + str(round(100*(ultimos >= 50).sum()/muestras,2)) + '%')

La probabilidad de ganar el juego es de 38.17%
