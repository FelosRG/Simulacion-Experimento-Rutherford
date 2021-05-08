# Simulacion Experimento Rutherford

Instrucciones:

Se requiere copilar primero el archivo Simulacion.f90 con f2py , f2py es un copilador especial de fortran que nos perimite crear una liberia con las subrutinas de Fortran. De esta forma accedemos desde python a la alta velocidad de computo de Fortran.

f2py -c Simulacion.f90 -m  ExperimentoRutherford

Una vez copilado el código ya es posible usar el jupyter notebook de manera usual.
