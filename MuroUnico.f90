!Programa de cálculo de un muro compuesto con materiales diversos, incluido PCM
!Pensado para cargar distintas configuraciones desde un archivo y analizarlas
!todas en lote
!Los cálculos de los muros se hacen con el método explícito
!Archivos usados en el programa:
!'configuraciones.txt' Archivo con las configuraciones a ensayar
!
!Variables usadas en el programa:
!ns,s: enteros con el número total de configuraciones, y el contador de configuraciones
!id: entero, identificador de configuración. Va a ser un código de 4 dígitos.
    !Se lee del archivo de configuraciones
!ncapa: número de capas del muro a cargar
!i,j:contador
!d,ho,t: contador de tiempo. d cuenta el día del año, ho, cuenta la hora del día, t el intervalo temporal
!hora: número de datos temporales incluidos en una hora
!contmin: contador de minutos, para grabar sólo los valores de tiempo cada minuto
!tipo(ncapa): código con el tipo de capa: 0 convección; 1 convencional; 2 PCM; 3 camara de aire
!nx(ncapa): número de nodos de capa capa del muro
!ic(ncapa-1): vector con la posición de las intercaras. Así aligero mucho el código de los cálculos de temperatura
    !En el elemento ic(i) aparece el número del nodo que corresponde a la intercara entre las capas i e i+1
!dt: intervalo de tiempo
!dx: intervalo espacial
!const: el valor dt/dx**2 se repite en muchos cálculos. Se calcula una vez, y no cambia
!tol: toleracia para el cálculo de la intercara convencional-PCM
!tintant,tintnuevo: valores para la temperatura interior, en el intante previo y en el actual
!densroom: densidad equivalente de la habitación
!densair: densidad del aire
!lroom: profundidad de la habitación
!cpair: calor específico del aire
!ren: renovaciones de aire
!muro: nombre del archivo con la configuración del muro
!kPCM: matriz con las propiedades de conductividad del PCM. 1a columna es ks, 2a columna Ts
    !y en 2a fila, de nuevo, primero kl y luego Tl
!long(ncapa),kond(ncapa),dens(ncapa),cp(ncapa),Rt(ncapa): vectores con las propiedades
    !termofísicas de las capas: espesor, conductividad, densidad, calor esp y resistencia de las capas de aire
!difus(ncapa): vector con las difusividades, para no tener que estar recalculandolas en cada iteración
!TempH(contmin,2): matriz para meter los datos de temperatura-entalpía
!knPCM(:): vector con los valores de conductividad nodal del instante anterior para la capa PCM
!h(:),Dh(:): vectores con el valor de entalpía y diferencia de entalpía nodales
!textair(hora),text(hora),qrad(hora): condiciones de contorno exteriores, temperatura sol-air y radiación efectiva en el interior
!qacond(hora): valor de la cantidad de energía de acondicionamiento. Se va a guardar la suma de todos los valores
    !de cada intervalo
!tempant(:),tempnuevo(:):vectores con valores de temperatura nodales del instante anterior y del nuevo
!condext:nombre del archivo con las condiciones exteriores
!capa: nombre del archivo con las propiedades de la capa
!filetemp,fileflux: nombre de los archivos de salida que voy a usar. Tienen que llevar el ID de la simulación

program MuroUnico
    use props
    implicit none

    integer::ns,s
    integer::id
    integer::ncapa
    integer::i,j,d,ho,t
    integer::hora
    integer::contmin
    integer,allocatable::tipo(:)
    integer,allocatable::nx(:)
    integer,allocatable::ic(:)

    real::dt,dx,const
    real::tol
    real::tintant,tintnuevo
    real::densroom,lroom,cpair,densair,ren
    !   real::kPCM(2,2) Valores asignados en el módulo props, como ks,Ts,kl y Tl
    real,allocatable::long(:),kond(:),dens(:),cp(:),Rt(:)
    real,allocatable::difus(:)
    !   real,allocatable::TempH(:,:)
    real,allocatable::knPCM(:),h(:),Dh(:)
    real,allocatable::textair(:),text(:),qrad(:)
    real,allocatable::qacond(:)
    real,allocatable::tempant(:),tempnuevo(:)

    character(len=20)::muro,condext
    character(len=20)::entalpia
    character(len=20)::capa
    character(len=14)::filetemp
    character(len=15)::fileflux

    !Datos iniciales constantes
    dt=0.01
    dx=0.005
    const=dt/dx**2
    tol=0.01
    hora=int(3600/dt+1)
    densroom=22.65
    lroom=3
    cpair=1007
    densair=1.2
    ren=1.5/3600

    !Carga de archivo de configuraciones e inicio de iteración en experimentos
    open(1,file='configuraciones.txt',status='old')
    read(1,*)ns !El archivo empieza con el número de configuraciones a simular
    do s=1,ns
        write(*,*)'Muro',s,'de',ns
        !Leer configuraciones
        read(1,*)id,muro,condext

        !Cargar los datos del muro
        open(2,file=muro,status='old')
        read(2,*)ncapa
        allocate(tipo(ncapa))
        allocate(long(ncapa),kond(ncapa),dens(ncapa),cp(ncapa),Rt(ncapa))
        allocate(difus(ncapa))
        long=0.
        kond=0.
        dens=0.
        cp=0.
        Rt=0.
        difus=0.
        do i=1,ncapa
            read(2,*)capa
            open(3,file=capa,status='old') !Abrir archivo con los datos de la capa
            read(3,*)tipo(i)
            select case (tipo(i))
                case(0)
                    read(3,*)Rt(i)
                case(1)
                    read(3,*)long(i)
                    read(3,*)kond(i)
                    read(3,*)dens(i)
                    read(3,*)cp(i)
                    difus(i)=kond(i)/dens(i)/cp(i)
                case(2)
                    read(3,*)long(i)
                    read(3,*)ks,Ts
                    read(3,*)kl,Tl
                    read(3,*)dens(i)
                    read(3,*)entalpia
                    !Carga de la matriz de temperatura-entalpía. Si se va a usar otro método de cálculo, se puede quitar este trozo
                    !La primera fila incluye el número de datos de entalpía a leer
                    open(4,file=entalpia,status='old')
                    read(4,*)contmin    !uso este contador porque todavía no lo estoy usando, por no generar otra varable
                    allocate(TempH(contmin,2))
                    do j=1,contmin
                        read(4,*)TempH(j,1),TempH(j,2)
                    end do
                    close(4)
                case(3)
                    read(3,*)long(i)
                    read(3,*)Rt(i)
            end select
            close(3)
        enddo
        close(2)

        !Cargar condiciones de contorno
        allocate(textair(hora),text(hora),qrad(hora))
        allocate (qacond(hora))
        open(2,file=condext,status='old')
        read(2,*)textair(1),text(1),qrad(1)
        textair(:)=textair(1)
        text(:)=text(1)
        qrad(:)=qrad(1)
        qacond=0.

        tintant=18.         !Inicializo las condiciones interiores a 18ºC

        !Inicializar condiciones iniciales. Para eso, crear primero vector de temperaturas
        allocate(nx(ncapa))
        forall (i=1:ncapa) nx(i)=int(long(i)/dx+1)
        allocate(ic(ncapa-1))
        forall (i=1:ncapa-1) ic(i)=sum(nx(1:i))-i+1
        allocate(tempant(ic(ncapa-1)),tempnuevo(ic(ncapa-1)))
        forall (i=1:ic(ncapa-1)) tempant(i)=text(1)-(text(1)-tintant)/(ic(ncapa-1)+1)*i
        tempnuevo=tempant
        allocate(knPCM(ic(ncapa-1)))
        allocate(h(ic(ncapa-1)))
        allocate(Dh(ic(ncapa-1)))
        knPCM=0.
        h=0.
        Dh=0.

        !Empezar algoritmo. Abrir archivo de salida de temperatura, y grabar instante inicial
        write(filetemp,'(A,I4.4,A)')'output',id,'.txt'
        write(fileflux,'(A,I4.4,A)')'outflux',id,'.txt'
        open(3,file=filetemp)
        write(3,*)'ID Ensayo',id !Escribo una cabecera con información
        do i=1,ncapa
            write(3,*)'Capa ',i, 'n nodos ',nx(i), 'Tipo',tipo(i)
        enddo
        write(3,*)1,0,text(1),tempant,tintant,qacond(1)
        contmin=1

        !Empezar iteración de tiempo
        do d=1,365
            do ho=0,23
                !Asigno como valor inicial el final de la hora anterior, leo el siguiente, e interpolo linealmente
                textair(1)=textair(hora)
                text(1)=text(hora)
                qrad(1)=qrad(hora)
                read(2,*)textair(hora),text(hora),qrad(hora)
                forall (i=2:hora-1)
                    textair(i)=textair(1)+(textair(hora)-textair(1))/3600*dt*(i-1)
                    text(i)=text(1)+(text(hora)-text(1))/3600*dt*(i-1)
                    qrad(i)=qrad(1)+(qrad(hora)-qrad(1))/3600*dt*(i-1)
                end forall
                qacond(:)=0.
                do t=2,hora
                    !Empezar a iterar por capas e intercapas
                    !Como paso cero, cargo los vectores de entalpía y knPCM para el material PCM
                    do i=1,ncapa
                        if (tipo(i)==2) then
                            do j=ic(i-1),ic(i) !busco los nodos entre las capas i-1 (ic(i-1)) e i+1 (ic(i))
                                knPCM(j)=k(tempant(j),tempant(j+1))
                                h(j)=HT(tempant(j))
                            enddo
                        endif
                    enddo

                    !Primero convección exterior
                    select case(tipo(2))
                        case(1)
                            tempnuevo(1)=tempant(1)+2/Rt(1)/dens(2)/cp(2)*dt/dx*(text(t-1)-tempant(1))+2*difus(2)*const*&
                                (tempant(2)-tempant(1))
                        case(2)
                            Dh(1)=2/Rt(1)/dens(2)*dt/dx*(text(t-1)-tempant(1))+2*knPCM(1)/dens(2)*const*(tempant(2)-tempant(1))
                            h(1)=h(1)+Dh(1)
                            tempnuevo(1)=TH(h(1))
                    end select

                    !Ahora nodos interiores, con un CASE.
                    do i=2,ncapa-1
                        select case(tipo(i))
                            case(1)
                                forall (j=ic(i-1)+1:ic(i)-1)
                                    tempnuevo(j)=tempant(j)+difus(i)*const*(tempant(j-1)-2*tempant(j)+tempant(j+1))
                                end forall
                            case(2)
                                forall (j=ic(i-1)+1:ic(i)-1)
                                    Dh(j)=const/dens(i)*(knPCM(j-1)*(tempant(j-1)-tempant(j))+knPCM(j)*(tempant(j+1)-tempant(j)))
                                end forall
                                h=h+Dh
                                do j=ic(i-1)+1,ic(i)-1
                                    tempnuevo(j)=TH(h(j))
                                enddo
                        end select
                    enddo

                    !Luego, el resto de intercapas, con un CASE para las combinaciones
                    do i=2,ncapa-2
                        select case(tipo(i))
                            case(1)     !convencional
                                select case(tipo(i+1))
                                    case(1) !contacto entre dos convencionales
                                        tempnuevo(ic(i))=tempant(ic(i))+2/(dens(i)*cp(i)+dens(i+1)*cp(i+1))*const*(kond(i)*&
                                            (tempant(ic(i)-1)-tempant(ic(i)))+kond(i+1)*(tempant(ic(i)+1)-tempant(ic(i))))
                                    case(2) !contacto convencional-PCM
                                        tempnuevo(ic(i))=tempant(ic(i))
                                        do
                                            Dh(ic(i))=2*const/dens(i+1)*(kond(i)*(tempant(ic(i)-1)-tempant(ic(i)))&
                                                +knPCM(ic(i))*(tempant(ic(i)+1)-tempant(ic(i)))&
                                                -dens(i)*cp(i)/2/const*(tempnuevo(ic(i))-tempant(ic(i))))
                                            h(ic(i))=h(ic(i))+Dh(ic(i))
                                            if(abs(tempnuevo(ic(i))-TH(h(ic(i))))<tol) then
                                                exit
                                            else
                                                tempnuevo(ic(i))=TH(h(ic(i)))
                                            endif
                                        enddo
                                        tempnuevo(ic(i))=TH(h(ic(i)))
                                    case(3) !contacto convencional-cámara de aire
                                        tempnuevo(ic(i))=tempant(ic(i))+2*difus(i)*const*(tempant(ic(i)-1)-&
                                            tempant(ic(i)))+2*dt/dens(i)/cp(i)/dx*(tempant(ic(i+1))-&
                                            tempant(ic(i)))/Rt(i+1)
                                end select
                            case(2)     !PCM
                                select case(tipo(i+1))
                                    case(1) !contacto PCM-convencional
                                        tempnuevo(sum(nx(1:i))-i+1)=tempant(sum(nx(1:i))-i+1)
                                        h(sum(nx(1:i))-i+1)=HT(tempant(sum(nx(1:i))-i+1))
                                        do
                                            Dh(ic(i))=2*const/dens(i)*(kond(i+1)*(tempant(ic(i)+1)-tempant(ic(i)))&
                                                +knPCM(ic(i)-1)*(tempant(ic(i)-1)-tempant(ic(i)))&
                                                -dens(i+1)*cp(i+1)/2/const*(tempnuevo(ic(i))-tempant(ic(i))))
                                            h(ic(i))=h(ic(i))+Dh(ic(i))
                                            tempnuevo(ic(i))=TH(h(ic(i)))
                                            if(abs(tempnuevo(ic(i))-TH(h(ic(i))))<tol) then
                                                exit
                                            else
                                                tempnuevo(ic(i))=TH(h(ic(i)))
                                            endif
                                        enddo
                                        tempnuevo(ic(i))=TH(h(ic(i)))
                                    !case(2)!no considero contacto entre dos PCMs
                                    case(3) !contacto PCM-cámara de aire
                                        Dh(ic(i))=2*knPCM(ic(i)-1)/dens(i)*const*(tempant(ic(i)-1)-&
                                            tempant(ic(i)))+2*dt/dens(i)/dx*(tempant(ic(i))-&
                                            tempant(ic(i+1)))/Rt(i+1)
                                        h(ic(i))=h(ic(i))+Dh(ic(i))
                                        tempnuevo(ic(i))=TH(h(ic(i)))
                                end select
                            case(3)     !contacto posterior de la cámara de aire
                                select case(tipo(i+1))
                                    case(1) !contacto cámara de aire-convencional
                                        tempnuevo(ic(i))=tempant(ic(i))+2*difus(i+1)*const*(tempant(ic(i)+1)-&
                                            tempant(ic(i)))+2*dt/dens(i+1)/cp(i+1)/dx*(tempant(ic(i-1))-&
                                            tempant(ic(i)))/Rt(i)
                                    case(2) !contacto cámara de aire-PCM
                                        Dh(ic(i))=2*knPCM(ic(i))/dens(i+1)*const*(tempant(ic(i)+1)-&
                                            tempant(ic(i)))+2*dt/dens(i+1)/dx*(tempant(ic(i-1))-&
                                            tempant(ic(i)))/Rt(i)
                                        h(ic(i))=h(ic(i))+Dh(ic(i))
                                        tempnuevo(ic(i))=TH(h(ic(i)))
                                    !case(3) no hay contacto entre dos cámaras de aire
                                end select
                        end select
                    enddo
                    !Convección interior
                    select case(tipo(ncapa-1))
                        case(1)
                            tempnuevo(ic(ncapa-1))=tempant(ic(ncapa-1))+2*difus(ncapa-1)*const*(tempant(ic(ncapa-1)-1)-&
                                tempant(ic(ncapa-1)))+2*dt/dx/cp(ncapa-1)/dens(ncapa-1)/Rt(ncapa)*(tintant-&
                                tempant(ic(ncapa-1)))
                        case(2)
                            Dh(ic(ncapa-1))=2/Rt(ncapa)/dens(ncapa-1)*dt/dx*(tintant-tempant(ic(ncapa-1)))+&
                                2*knPCM(ic(ncapa-1)-1)/dens(ncapa-1)*const*(tempant(ic(ncapa-1)-1)-tempant(ic(ncapa-1)))
                            h(ic(ncapa-1))=h(ic(ncapa-1))+Dh(ic(ncapa-1))
                            tempnuevo(ic(ncapa-1))=TH(h(ic(ncapa-1)))
                    end select

                    !Balance de aire interior
                    tintnuevo=tintant+dt/densroom/lroom/cpair*((tempnuevo(ic(ncapa-1))-tintant)/Rt(ncapa)+qrad(t)-&
                    lroom*ren*densair*cpair*(tintant-textair(t-1)))
                    if ((tintnuevo<20).and.((d<=130).or.(d>305))) then
                        if(((ho.lt.7).or.(ho==23)).and.tintnuevo<17) then
                            tintnuevo=17
                            qacond(t)=densroom*lroom*cpair*(tintant-tintnuevo)/dt-(tempnuevo(ic(ncapa-1))-tintant)/Rt(ncapa)-&
                            qrad(t)+lroom*ren*densair*cpair*(tintant-textair(t-1))
                        elseif((ho.ge.7).and.(ho<23)) then
                            tintnuevo=20
                            qacond(t)=densroom*lroom*(tintant-tintnuevo)/dt-(tempnuevo(ic(ncapa-1))-tintant)/Rt(ncapa)-qrad(t)+&
                            lroom*ren*densair*cpair*(tintant-textair(t-1))
                        endif
                    endif
                    if ((tintnuevo>25).and.((d>130).and.(d<=305))) then
                        tintnuevo=25
                        qacond(t)=densroom*lroom*cpair*(tintant-tintnuevo)/dt-(tempnuevo(ic(ncapa-1))-tintant)/Rt(ncapa)-qrad(t)+&
                        lroom*ren*densair*cpair*(tintant-textair(t-1))
                    endif

                    !Calculados todos los valores de temperatura, se guardan en el vector tempant para la siguiente iteración
                    tempant=tempnuevo
                    tintant=tintnuevo
                    !Comprobación de grabado en archivo de las temperaturas y cambio de vectores nuevos al antiguo. Cada 15 mins
                    if(contmin==900/dt) then    !Si cambio el intervalo, tendré que cambiar el cálculo de flujos
                        write(3,*)d,ho,text(t),tempnuevo,tintnuevo,qacond(t)
                        contmin=0
                    else
                        contmin=contmin+1
                    endif
                enddo
            enddo
        enddo !fin de la iteración en el tiempo
        close(3)
        close(2)

        !Cálculo de los flujos de calor en las superficies. Cierre de archivos y deallocates de funciones
        !Necesito un vector especial para las temperaturas que voy leyendo, y otro para sacar las temperaturas
        !de los extremos sobre las que calcular los flujos de calor
        !En lugar de crear variables nuevas, uso variables ya definidas: tempant para la lectura de datos,
        !tempnuevo para los valores de los extremos, y text para los flujos de calor calculados
        deallocate(textair,text,qrad,qacond)
        deallocate(tempant,tempnuevo)

        allocate(tempant(ic(ncapa-1)+5))  !Además de los nodos del material hay que leer también las temps de contorno
        allocate(tempnuevo(4))  !las dos temperaturas exteriores y las dos superficiales
        allocate(text(2))       !Los flujos, sólo son dos, cada uno en una cara

        open(2,file=filetemp)
        do t=1,ncapa+1
            read(2,*)
        enddo
        open(3,file=fileflux)

        do t=1,int(365*24*4) !El último multiplicando, el 4, hace referencia a que estoy grabando datos cada cuarto
            read(2,*)tempant    !de hora. Si cambio el grabado, lo tendré que cambiar
            tempnuevo(1)=tempant(3)
            tempnuevo(2)=tempant(4)
            tempnuevo(3)=tempant(ic(ncapa-1)+3)
            tempnuevo(4)=tempant(ic(ncapa-1)+4)
            text(1)=(tempnuevo(1)-tempnuevo(2))/Rt(1) !positivo del exterior al interior
            text(2)=(tempnuevo(3)-tempnuevo(4))/Rt(ncapa) !positivo del exterior al interior
            write(3,*)text
        enddo
        close(2)
        close(3)

        deallocate(tipo)    !desalojo el resto de variables porque en la siguiente configuración pueden cambiar dimensiones
        deallocate(long,kond,dens,cp,Rt)
        deallocate(difus)
        if (allocated(TempH)) deallocate(TempH)
        deallocate(nx)
        deallocate(ic)
        deallocate(knPCM)
        deallocate(h,Dh)
        deallocate(text,tempant,tempnuevo)

    enddo !fin de la iteración en configuraciones
    close(1)
end program MuroUnico
