
! **************************************************************************!
! IMPORTANTE: NO QUITAR LOS COMENTARIOS "!f2py intent(in,out) :: ... "      !
!             en los inicios de la linea, es una orden especial dada al     !
!             copilador de Fortran para Python.                             !
! **************************************************************************

subroutine distancia2(V1,V2,distancia)
    ! Calcula la distancia al cuadrado de dos puntos.

    ! Posición de la partícula      (INPUT)
    real*8 , intent(in)  :: V1(0:2) , V2(0:2)
    ! Distancia entre los vectores  (OUTPUT)
    real*8 , intent(out) :: distancia
    distancia =  sum((V1 - V2)**2)
end subroutine

subroutine NearestAtomCP(X,MAX_DESLOCACION,DIAMETRO_CELDA,POS_ATOMO)
    ! Encuentra al núcleo más cercano dentro de la red cristalina.

    real*8 , intent(in)   ::          X(0:2)             ! Posición de la partícula alfa (Input) 
    real*8 , intent(out)  ::  POS_ATOMO(0:2)             ! Posición del nucleo           (Output)
    ! Parámetros de configuración
    real*8 , intent(in)   :: DIAMETRO_CELDA              ! Diámetro de la celda unitaria.
    real*8 , intent(in)   :: MAX_DESLOCACION             ! Máxima fracción del diámetro de la celda en la que puede desplazarse cada núcleo. 
    ! Otras Variables...
    real*8  :: deslocacion
    ! Vectores       
    real*8  :: X_NORM(0:2)             ! Se noramliza la posición respecto al diámetro de la celda unitaria.
    real*8  :: POS_CENTRO_CELDA(0:2)
    real*8  :: NUM_CELDA(0:2)
    real*8  :: SIGNO(0:2)              ! Matriz que dererminará

    ! Paso 1:  Obtenemos la celda del grid en la que está la partícula alfa.
    X_NORM    = X / DIAMETRO_CELDA
    NUM_CELDA = nint(X_NORM)
    ! Paso 2: Obtenemos la posición del centro de la celda en la que se encuentra la partícula alfa.
    POS_CENTRO_CELDA  = NUM_CELDA*DIAMETRO_CELDA
    ! Paso 3: Mediante condiciones lógicas, obtenemos los 4 átomos más cercanos al átomo dentro de la celda unitaria.
    ! Virtualmente dividimos el cubo de la celda unitaria en 8 subregiones, cada subregión tiene asociado 4 átomos del cristal.
    do i=0,2
        ! Asignación de la matriz de condición para cada componente.
        if (X(i) > POS_CENTRO_CELDA(i)) then
            SIGNO(i) = 1
        else
            SIGNO(i) = -1
        end if
    end do
    ! Este es el átomo más cercano.
    POS_ATOMO =   POS_CENTRO_CELDA + SIGNO*(DIAMETRO_CELDA/2)
    ! Agregamos la deslocación aleatoria.
    call random_number(deslocacion)
    deslocacion = (deslocacion*2 - 1)*(DIAMETRO_CELDA*MAX_DESLOCACION/2)  ! Cambiamos el dominio de los números aleatorios, de 0-1 a -1,1. Y escalamos.
    POS_ATOMO = POS_ATOMO + deslocacion
    ! Definimos que la placa empieza en el 0 de la coordenada Y, por lo que nos aseguramos de que todos los átomos esten en coordenadas postivas.
    POS_ATOMO(1) = abs(POS_ATOMO(1))
end subroutine

subroutine NearestAtomBCC(X,MAX_DESLOCACION,DIAMETRO_CELDA,POS_ATOMO)
    ! Encuentra al núcleo más cercano dentro de la red cristalina.

    real*8 , intent(in)   ::          X(0:2)             ! Posición de la partícula alfa (Input) 
    real*8 , intent(out)  ::  POS_ATOMO(0:2)             ! Posición del nucleo           (Output)
    ! Parámetros de configuración
    real*8 , intent(in)   :: DIAMETRO_CELDA              ! Diámetro de la celda unitaria.
    real*8 , intent(in)   :: MAX_DESLOCACION             ! Máxima fracción del diámetro de la celda en la que puede desplazarse cada núcleo. 
    ! Otras Variables...
    real*8  :: distancia , distancia2  ! Variables que se relacionan con distancias.
    ! Vectores          
    real*8  :: X_NORM(0:2)             ! Se noramliza la posición respecto al diámetro de la celda unitaria.
    real*8  :: POS_CENTRO_CELDA(0:2)
    real*8  :: NUM_CELDA(0:2)
    real*8  :: SIGNO(0:2)              ! Matriz que dererminará
    ! Matrices
    real*8  :: deslocacion(0:1,0:2)
    real*8  :: celdas(0:1,0:2)

    ! Paso 1:  Obtenemos la celda del grid en la que está la partícula alfa.
    X_NORM    = X / DIAMETRO_CELDA
    NUM_CELDA = nint(X_NORM)
    ! Paso 2: Obtenemos la posición del centro de la celda en la que se encuentra la partícula alfa.
    POS_CENTRO_CELDA  = NUM_CELDA*DIAMETRO_CELDA
    ! Paso 3: Mediante condiciones lógicas, obtenemos los 4 átomos más cercanos al átomo dentro de la celda unitaria.
    ! Virtualmente dividimos el cubo de la celda unitaria en 8 subregiones, cada subregión tiene asociado 4 átomos del cristal.
    do i=0,2
        ! Asignación de la matriz de condición para cada componente.
        if (X(i) > POS_CENTRO_CELDA(i)) then
            SIGNO(i) = 1
        else
            SIGNO(i) = -1
        end if
    end do
    ! Agregamos el primer átomo de los 4 a la lista. Este es el átomo de la esquina de lla subregión.
    celdas(0,:) =  POS_CENTRO_CELDA + SIGNO*(DIAMETRO_CELDA/2)
    celdas(1,:) =  POS_CENTRO_CELDA
    ! Agregamos la deslocación aleatoria.
    call random_number(deslocacion)
    deslocacion = (deslocacion*2 - 1)*(DIAMETRO_CELDA*MAX_DESLOCACION/2)  ! Cambiamos el dominio de los números aleatorios, de 0-1 a -1,1. Y escalamos.
    celdas = celdas + deslocacion
    ! Encontramos el átomo más cercano.
    distancia = sum((celdas(0,:) - X)**2)
    POS_ATOMO = celdas(0,:)
    distancia2 = sum((celdas(1,:) - X)**2)
    if (distancia2 < distancia) then
        distancia = distancia2
        POS_ATOMO = celdas(1,:)
    end if
    ! Definimos que la placa empieza en el 0 de la coordenada Y, por lo que nos aseguramos de que todos los átomos esten en coordenadas postivas.
    POS_ATOMO(1) = abs(POS_ATOMO(1))
end subroutine

subroutine NearestAtomFCC(X,MAX_DESLOCACION,DIAMETRO_CELDA,POS_ATOMO)
    ! Encuentra al núcleo más cercano dentro de la red cristalina.

    real*8 , intent(in)   ::          X(0:2)             ! Posición de la partícula alfa (Input) 
    real*8 , intent(out)  ::  POS_ATOMO(0:2)             ! Posición del nucleo           (Output)
    ! Parámetros de configuración
    real*8 , intent(in)   :: DIAMETRO_CELDA              ! Diámetro de la celda unitaria.
    real*8 , intent(in)   :: MAX_DESLOCACION             ! Máxima fracción del diámetro de la celda en la que puede desplazarse cada núcleo. 
    ! Otras Variables...
    real*8  :: distancia , distancia2  ! Variables que se relacionan con distancias.
    ! Vectores
    real*8  :: BASE(0:2)           
    real*8  :: X_NORM(0:2)             ! Se noramliza la posición respecto al diámetro de la celda unitaria.
    real*8  :: POS_CENTRO_CELDA(0:2)
    real*8  :: NUM_CELDA(0:2)
    real*8  :: SIGNO(0:2)              ! Matriz que dererminará
    ! Matrices
    real*8  :: deslocacion(0:3,0:2)
    real*8  :: celdas(0:3,0:2)

    ! Paso 1:  Obtenemos la celda del grid en la que está la partícula alfa.
    X_NORM    = X / DIAMETRO_CELDA
    NUM_CELDA = nint(X_NORM)
    ! Paso 2: Obtenemos la posición del centro de la celda en la que se encuentra la partícula alfa.
    POS_CENTRO_CELDA  = NUM_CELDA*DIAMETRO_CELDA
    ! Paso 3: Mediante condiciones lógicas, obtenemos los 4 átomos más cercanos al átomo dentro de la celda unitaria.
    ! Virtualmente dividimos el cubo de la celda unitaria en 8 subregiones, cada subregión tiene asociado 4 átomos del cristal.
    do i=0,2
        ! Asignación de la matriz de condición para cada componente.
        if (X(i) > POS_CENTRO_CELDA(i)) then
            SIGNO(i) = 1
        else
            SIGNO(i) = -1
        end if
    end do
    ! Agregamos el primer átomo de los 4 a la lista. Este es el átomo de la esquina de lla subregión.
    celdas(0,:) =   POS_CENTRO_CELDA + SIGNO*(DIAMETRO_CELDA/2)
    ! Agregamos los otros tres átomos.
    do i=0,2
        BASE(0:2)     = (/0.0,0.0,0.0/)
        BASE(i)       = SIGNO(i)*DIAMETRO_CELDA/2
        celdas(i+1,:) = POS_CENTRO_CELDA + BASE
    end do
    ! Agregamos la deslocación aleatoria.
    call random_number(deslocacion)
    deslocacion = (deslocacion*2 - 1)*(DIAMETRO_CELDA*MAX_DESLOCACION/2)  ! Cambiamos el dominio de los números aleatorios, de 0-1 a -1,1. Y escalamos.
    celdas = celdas + deslocacion
    ! Encontramos el átomo más cercano.
    distancia = sum((celdas(0,:) - X)**2)
    POS_ATOMO = celdas(0,:)
    do i=1,3
        distancia2 = sum((celdas(i,:) - X)**2)
        if (distancia2 < distancia) then
            distancia = distancia2
            POS_ATOMO = celdas(i,:)
        end if
    end do
    ! Definimos que la placa empieza en el 0 de la coordenada Y, por lo que nos aseguramos de que todos los átomos esten en coordenadas postivas.
    POS_ATOMO(1) = abs(POS_ATOMO(1))
end subroutine

subroutine NearestAtom(CATEGORIA,X,MAX_DESLOCACION,DIAMETRO_CELDA,POS_ATOMO)
    ! Realiza la selección de la estructura cristalina.

    integer , intent(in)  :: CATEGORIA
    real*8 , intent(in)   ::          X(0:2)             ! Posición de la partícula alfa (Input) 
    real*8 , intent(out)  ::  POS_ATOMO(0:2)             ! Posición del nucleo           (Output)
    ! Parámetros de configuración
    real*8 , intent(in)   :: DIAMETRO_CELDA              ! Diámetro de la celda unitaria.
    real*8 , intent(in)   :: MAX_DESLOCACION             ! Máxima fracción del diámetro de la celda en la que puede desplazarse cada núcleo. 
    
    if    (CATEGORIA == 1) then
       call NearestAtomCP(X,MAX_DESLOCACION,DIAMETRO_CELDA,POS_ATOMO)
    elseif(CATEGORIA == 2) then
        call NearestAtomBCC(X,MAX_DESLOCACION,DIAMETRO_CELDA,POS_ATOMO)
    elseif(CATEGORIA == 3) then
        call NearestAtomFCC(X,MAX_DESLOCACION,DIAMETRO_CELDA,POS_ATOMO)
    ! Si no se entrega un modo valido, se asume una estructura primitiva cúbica.
    else
        call NearestAtomCP(X,MAX_DESLOCACION,DIAMETRO_CELDA,POS_ATOMO)
    end if
end subroutine

subroutine EspacioInteratomico(X,V,ATOM_RADIUS)
    ! Avanza a la partícula cuando está en el espacio interátomico.
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15)
    ! Posición y velocidad de entrada - salida. (INOUT)
    real*8 , intent(inout) :: X(0:2)
!f2py intent(in,out) :: X
    real*8 , intent(in)    :: V(0:2)
    real*8 , intent(in)    :: ATOM_RADIUS
    real*8 :: dt
    real*8 :: norma_velocidad
    real*8 :: V_Zero(0:2) = (/0.0_dp , 0.0_dp , 0.0_dp/)
    ! Paso 1: Encontramos la magnitud de la velocidad.
    call distancia2(V,V_Zero,norma_velocidad)
    dt = 0.25*ATOM_RADIUS / sqrt(norma_velocidad) ! Avanza solo una fracción del radio atómico.
    ! Paso 2: Actualizamos el valor de la posición.
    X = X + dt*V
end subroutine

subroutine MetodoNumerico(X,XC,V,Z_A,ATOM_RADIUS,metrica_num_pasos,metrica_num_interacciones,metrica_desviacion,VEL_INICIAL,norma_v)
    ! Subroutina que se encarga en su totalidad de ejecutar los métodos numéricos.

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    ! Valores de entrada y salida.
    real*8 , intent(in) :: XC(0:2)
    integer, intent(in) :: Z_A
    real*8 , intent(in) :: ATOM_RADIUS
    real*8 , intent(in) :: VEL_INICIAL
    real*8 , intent(inout) :: X(0:2)
!f2py intent(in,out) :: X
    real*8 , intent(inout) :: V(0:2)
!f2py intent(in,out) :: V
    real*8 , intent(out) :: norma_v
    ! Constantes
    real*8 , parameter  :: ALPHA_MASS   = 6.645e-27_dp!_dpr         ! [ kg ]
    real*8 , parameter  :: UNIT_CHARGE  = 1.6021e-19_dp!_dpr        !-19 [ C ]
    real*8 , parameter  :: K_CONSTANT   = 8.998e09_dp!_dpr         ! [ N m^2 / C ]
    ! Variables como métricas.
    logical          , intent(inout) :: metrica_desviacion
    integer          , intent(inout) :: metrica_num_pasos
    integer          , intent(inout) :: metrica_num_interacciones
    ! Otras variables.
    real*8  :: Radio_interaccion , Radio_interaccion2
    real*8  :: Carga_atomo , Carga_alpha
    real*8  :: dt , distancia
    real*8  :: var_a,var_phi
    real*8  :: multiplicador
    ! Más vectores.
    real*8  :: A(0:2)
    real*8  :: V_ZERO(0:2) = (/0.0_dp , 0.0_dp , 0.0_dp/)

    Radio_interaccion  = 0.75*ATOM_RADIUS
    Radio_interaccion2 = Radio_interaccion**2
    call distancia2(X,XC,distancia)

    ! EXPERIMENTAL : ABSORCION DE ENERGÍA.
    call distancia2(V,V_ZERO,norma_v)
    norma_v = sqrt(norma_v)
    multiplicador = VEL_INICIAL / norma_v

    if (distancia >=  Radio_interaccion2) then
        call EspacioInteratomico(X,V,ATOM_RADIUS)
        metrica_num_pasos = metrica_num_pasos + 1
    else
        metrica_desviacion = .TRUE.
        ! Cuerpo principal del programa.
        metrica_num_interacciones = metrica_num_interacciones + 1
        Carga_atomo = Z_A*UNIT_CHARGE
        Carga_alpha =  2.*UNIT_CHARGE 
        var_a          = K_CONSTANT*Carga_atomo*Carga_alpha /  ALPHA_MASS
        distancia = sqrt(distancia)
        do while(distancia < Radio_interaccion)
            ! Método de verlet con velocidades explicitas.
            dt = distancia*2.5e-09_dp*multiplicador ! El paso temporal es proporcional a la distancia del atómo con la que está interactuando.
            var_phi = var_a / distancia**3
            A = var_phi*(X-XC)
            ! Actualizamos posición
            X = X + dt*V
            call distancia2(X,XC,distancia) ! Actualización de la distancia.
            distancia = sqrt(distancia)
            var_phi = var_a / distancia**3
            ! Actualizamos velocidad.
            V  = V  + 0.5*dt*(A + var_phi*(X-XC))
            metrica_num_pasos =  metrica_num_pasos + 1
        end do
        V = V*0.9999
        ! Termina la subrutina cuando la partícula alfa sale del radio de interacción con el átomo.
    end if
end subroutine 

subroutine Simulacion(N,CATEGORIA,Z_A,NUM_LAMINAS,ATOM_RADIUS,DIAMETRO_CELDA,&
                      VEL_INICIAL,MAX_DESLOCACION,metrica_num_desviadas,&
                      metrica_num_detenidas,metrica_num_rebotadas,FILE_NAME,verbose)
    ! Ejecuta Toda la simulación de una sola partícula.

    INTEGER , PARAMETER  :: dp = SELECTED_REAL_KIND(15)
    ! INPUTS
    integer , intent(in) :: N
    integer , intent(in) :: CATEGORIA
    integer , intent(in) :: Z_A
    real*8  , intent(in) :: NUM_LAMINAS
    real*8  , intent(in) :: ATOM_RADIUS
    real*8  , intent(in) :: DIAMETRO_CELDA
    real*8  , intent(in) :: VEL_INICIAL
    real*8  , intent(in) :: MAX_DESLOCACION
    ! Otros Inputs
    logical , intent(in) :: verbose
    Character(len=40) , intent(in) :: FILE_NAME
    ! Matrices
    real*8 :: random(N,0:2)
    ! Vectores
    real*8  :: X(0:2)
    real*8  :: POS_ATOMO(0:2)
    real*8  :: V(0:2)
    real*8  :: V_ZERO(0:2) = (/0.0_dp , 0.0_dp , 0.0_dp/)
    ! Variables
    real*8, parameter    :: DIAMETRO_RAYO  = 100e-10 ! Ancho del rayo de la fuente radioactiva. (Nos asegura variabilidad.)
    integer              :: Atomos_Lamina  = 200
    ! Metricas
    integer :: metrica_num_pasos
    integer :: metrica_num_interacciones
    integer , intent(out) :: metrica_num_detenidas
    integer , intent(out) :: metrica_num_rebotadas
    integer , intent(out) :: metrica_num_desviadas
    ! Variables de estado (indica si pasó algun evento especial con la partícula)
    logical :: estado_detenido
    logical :: metrica_desviacion
    ! Otras variables
    real*8  :: ancho_capas
    real*8  :: norma_v
    real*8  :: norma
    integer :: porcentaje
    ! iniciamos las metricas
    metrica_num_detenidas = 0
    metrica_num_rebotadas = 0
    metrica_num_desviadas = 0
    ! Inicio de la simulación
    if (verbose)  print*, "Iniciando simulación ... "
    ! Paso 1: Calculamos variables, iniciamos archivo de datos.
    open(1,file=FILE_NAME)
    call random_number(random)
    random = random*2 -1 ! Cambiamos dominio a (-1 , 1)
    porcentaje = N / 10
    ancho_capas = Atomos_Lamina*DIAMETRO_CELDA*NUM_LAMINAS
    ! Paso 2: Bucle para la simulación de cada partícula.
    do i=1,N
        ! Reseteamos los parámetros de estado.
        estado_detenido    = .False.
        metrica_desviacion = .FALSE.
        ! Inicializamos los parámetros iniciales
        X(0:2) = (/0.0_dp , 0.0_dp , 0.0_dp/)
        V(0:2) = (/0.0_dp , VEL_INICIAL , 0.0_dp /)
        ! Obtenemos la posición aleatoria
        X    = random(i,:)*(DIAMETRO_CELDA + 1/1e12)
        X(1) = -0.9*ATOM_RADIUS
        ! Ciclo principal , se repite hasta salir de la placa.
        do while(abs(X(1)) < ancho_capas)
            !Paso 1 : Encontramos el atomo más cercano de la red cristalina.
            call NearestAtom(CATEGORIA,X,MAX_DESLOCACION,DIAMETRO_CELDA,POS_ATOMO)
            !Paso 2 : Interacción partícula alpha atomo
            call MetodoNumerico(X,POS_ATOMO,V,Z_A,ATOM_RADIUS,metrica_num_pasos,metrica_num_interacciones,&
            metrica_desviacion,VEL_INICIAL,norma_v)

            ! Acciones caundo la partícula sea detenida.
            if (norma_v < 0.05*VEL_INICIAL) then
                metrica_num_detenidas = metrica_num_detenidas + 1
                estado_detenido = .True.
                EXIT ! Rompemos el ciclo while.
            end if

        end do

        ! Saltamos ciclo si la partícula es detenida.
        !if (estado_detenido .eqv. .True.) CYCLE

        ! Contamos el número de rebote.
        if (V(1) < 0) metrica_num_rebotadas = metrica_num_rebotadas + 1

        ! Normalizamos el vector de velocidad
        call distancia2(V,V_ZERO,norma)
        norma = sqrt(norma)
        V     = V / norma
        ! Escribimos resultado. (Solo escribe el resultado de las partículas que fueron desviadas)
        if (metrica_desviacion .eqv. .TRUE.) then 
            write(1,*) V(0) , V(1) , V(2) , norma / VEL_INICIAL
            metrica_num_desviadas = metrica_num_desviadas + 1
        end if
        if (verbose) then
            if (mod(i,porcentaje)  == 0) then
                print*, (i*100) / N , "%"
            end if
        end if
    end do

    if (verbose) then
        print*, " Num de partículas simuladas : " , N
        print*, " Num de partículas desviadas : " , metrica_num_desviadas
        print*, " "
        print*, " Nucleos prom Interactuados  : " , metrica_num_interacciones / metrica_num_desviadas
        print*, " Pasos   prom efectuados     : " , metrica_num_pasos         / N
    end if 
endsubroutine


subroutine Prueba()
    real*8  :: a = 2
    real*16  :: b
    b = sqrt(a)
    print*, b
end subroutine


program main
    integer :: N  = 10000
    integer :: CATEGORIA = 3
    integer :: Z_A             = 79
    real*8  :: NUM_LAMINAS     = 10
    real*8  :: ATOM_RADIUS     = 1.44e-10
    real*8  :: DIAMETRO_CELDA  = 4.07e-10
    real*8  :: VEL_INICIAL     = 1.57e7
    real*8  :: MAX_DESLOCACION = 0.5
    integer :: metrica_num_desviadas
    integer :: metrica_num_detenidas
    integer :: metrica_num_rebotadas
    Character(len=40) :: FILE_NAME = "Prueba.txt"
    logical :: verbose = .TRUE.
    call Simulacion(N,CATEGORIA,Z_A,NUM_LAMINAS,ATOM_RADIUS,DIAMETRO_CELDA,VEL_INICIAL,&
                    MAX_DESLOCACION,metrica_num_desviadas,metrica_num_detenidas,metrica_num_rebotadas,FILE_NAME,verbose)
    print*, "Paradas  :",metrica_num_detenidas
    print*, "rebotadas:",metrica_num_rebotadas
end program 