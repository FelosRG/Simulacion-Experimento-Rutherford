#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  8 15:16:26 2021

@author: felos
"""
import ExperimentoRutherfordFortran
import numpy as np
import multiprocessing
import math

def SimulacionOneCPU(N        ,
               categoria      ,
               z_a            ,
               num_laminas    ,
               atom_radius    ,
               diametro_celda ,
               vel_inicial    , 
               max_deslocacion,
               clave_proceso  ):
    
    metrica_desviaciones = ExperimentoRutherfordFortran.simulacion(N,
                                                            categoria      ,
                                                            z_a            ,
                                                            num_laminas    ,
                                                            atom_radius    ,
                                                            diametro_celda ,
                                                            vel_inicial    ,
                                                            max_deslocacion,
                                                            f"Resultado{clave_proceso}.txt", 
                                                            verbose=False)
    return metrica_desviaciones
    
    
    


def ExperimentoNumerico(N,
                        cpus         = 2  ,
                        tipo_cristal ="fcc",
                        z_a          = 79 ,
                        num_laminas  = 1  ,
                        atom_radius  = 1.44e-10,
                        diametro_celda = 4.07e-10,
                        vel_inicial    = 1.57e7,
                        max_deslocacion= 0.5):
    
    
    # Hacemos parámetros compatibles.
    if tipo_cristal == "pc":
        categoria = 1
    elif tipo_cristal == "bcc":
        categoria = 2
    elif tipo_cristal == "fcc":
        categoria = 3
    else:
        raise ValueError("Los tipos de cristal disponibles son \"pc\",\"bcc\" y \"fcc\"")
    metrica = 0
    if cpus == 1:
        # Corremos en una sola cpu.
        metrica_desviaciones = SimulacionOneCPU(N              ,
                                                categoria      ,
                                                z_a            ,
                                                num_laminas    ,
                                                atom_radius    ,
                                                diametro_celda ,
                                                vel_inicial    , 
                                                max_deslocacion,
                                                0 )
        
        # Recopilamos Resultado
        array_resultado = np.loadtxt(f"Resultado{0}.txt")
        # Diccionario salida 
        Output = {"datos":array_resultado,"metrica":metrica}
        return Output
    else:
        # Repartimos el trabajo entre las cpus.
        tamaño_rebanada            = math.floor(N/cpus)
        processes = []
        for i in range(cpus):
            process = multiprocessing.Process(target=SimulacionOneCPU,args=(tamaño_rebanada,
                                                                            categoria      ,
                                                                            z_a            ,
                                                                            num_laminas    ,
                                                                            atom_radius    ,
                                                                            diametro_celda ,
                                                                            vel_inicial    ,
                                                                            max_deslocacion,
                                                                            i              ))
            processes.append(process)
            process.start()
        for process in processes:
            process.join()
            
        # Recopilamos los resultados.
        array_resultado = np.loadtxt(f"Resultado{0}.txt")
        for num_resultado in range(i):
            i += 1
            resultado       =  np.loadtxt(f"Resultado{num_resultado}.txt")
            try:
                array_resultado =  np.append(array_resultado,resultado,axis=0)
            except ValueError:
                pass
        
        # Diccionario de salida.
        # ! Falta sumar las desviaciones de las rebanadas.
        Output = {"datos":array_resultado,"metrica":metrica}
        return Output