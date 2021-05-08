# Simulacion Experimento Rutherford
![alt text](https://github.com/FelosRG/Simulacion-Experimento-Rutherford/blob/main/Figuras/Portada.png?raw=true)

## Instrucciones:

Se requiere copilar primero el archivo Simulacion.f90 con f2py , f2py es un copilador especial de fortran incluido en la libería de numpy que nos perimite crear una libería con las subrutinas de Fortran que se importan como funcione regulares en Python. 
De esta forma accedemos desde python a la alta velocidad de computo de Fortran.

### Ubuntu
Nos aseguramos primero de que tengamos numpy.
```
sudo apt-get install python-numpy
```
Copilamos el código con f2py
```
f2py -c Simulacion.f90 -m  ExperimentoRutherford
```

Una vez copilado el código ya es posible usar el jupyter notebook de manera usual.
