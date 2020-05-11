C      PROGRAMA cincotta1.f
 
C---------------------------------------------------------------------
C     SEGUN EL PAPER DE: Cincotta, M\'endez & N\u'\~nez 1995, ApJ 449, 231.
C
C
C pinicial  ES EL PERIODO INICIAL PARA LA BUSQUEDA.
C pfinal    ES EL PERIODO FINAL PARA LA BUSQUEDA. 
C nperiod   ES EL NUMERO DE PERIODOS QUE SE VAN A ITERAR.
C fa1 Y fa2 SON EL NUMERO DE BINES DE LA CURVA DE LUZ (EN LAS          C    DIRECCIONES 
C           phi Y u RESPECTIVAMENTE.
C---------------------------------------------------------------------
C      SI EL NUMERO DE DATOS DE LA CURVA DE LUZ ES MAYOR QUE 1000
C      EL TAMAN~O DE LAS TRES SIGUIENTES MATRICES DEBE SER CAMBIADO.
C---------------------------------------------------------------------
       IMPLICIT REAL*8(a-h,o-z)
       dimension t(1000),v(1000),f(1000)
       open(50,file='datos.dat',status='old')
       open(21,file='ceph.17.I .dat',status='old')
       open(51,file='ceph.17.IS.dat',status='unknown')

C      LECTURA DEL ARCHIVO DE PARAMETROS datos.dat
       read(50,*) pinicial 
       read(50,*) pfinal   
       read(50,*) nperiod        
       read(50,*) fa1      
       read(50,*) fa2      
       read(50,*) NDATOS   
       close(50)

C      CONTROL DEL NUMERO DE BINES       
C      EL NUMERO DE BINES NO DEBE EXCEDER (1000*1000) 
       if(fa1.gt.1000.or.fa2.gt.1000) then
         write(*,*) 'Number of bins should be < 1000!'
         write(*,*) 'Enter new number of bins (PHI, U):'
         read(*,*) fa1,fa2
         if(fa1.gt.1000.or.fa2.gt.1000) stop
       endif

C      dperiod   ES EL PASO DE PERIODO (UNIFORME).             
       dperiod=(pfinal-pinicial)/nperiod

C      BUSQUEDA DE LOS MINIMOS Y MAXIMOS VALORES DE LA CURVA DE LUZ.
       n  = 0
       tmin =  10000000.
       tmax = -10000000.
       vmin =  1000000.
       vmax = -1000000.

       DO 10 K=1,NDATOS 
       n=n+1
       read(21,*) t(n),v(n)
       tmin=min(tmin,t(n))
       tmax=max(tmax,t(n))
       vmin=min(vmin,v(n))
       vmax=max(vmax,v(n))
10     CONTINUE  
       ndata=n-1
       close(21)

C      ELECCION DEL TIEMPO 0 PARA CALCULAR LA FASE
       to = tmin
C      NORMALIZACION DE LA CURVA DE LUZ
       do n=1,ndata
         v(n)=(v(n)-vmin)/(vmax-vmin)
       enddo

C      ITERACION DEL  PERIODO ENTRE pinicial Y  pfinal      
       do k=0,nperiod
         p  = pinicial + (dfloat(k) * dperiod)
C      CALCULO DE LA FASE         
         do i=1,ndata
            f(i)=((t(i)-to)/p) - int((t(i)-to)/p)
C      CORRIMIENTO DE LA FASE SI ES < 0.            
            if(f(i).lt.0) f(i)=f(i)+1.0d0
         enddo
C      LLAMADO DE LA SUBRUTINA QUE CALCULA LA ENTROPIA         
         call entrop(f,v,s,ndata,fa1,fa2)
C      ESCRITURA DEL ARCHIVO DE SALIDA CON: PERIODO, ENTROPIA         
         write(51,151)p,s
 151     format(e14.8,1x,e14.8)         
       enddo
       close(51)
       stop
       end

c***********************************************************************
c SUBRUTINA QUE EVALUA LA ENTROPIA DE LA CURVA DE LUZ
c***********************************************************************
       SUBROUTINE entrop(x,y,s,ndata,fa1,fa2) 
       IMPLICIT REAL*8(A-H,O-Z)
       dimension his(1000,1000),x(*),y(*)
C      LA MATRIZ his(i,j) ALMACENA EL NUMERO DE DATOS EN EL BIN (i,j)

C      SE REALIZA LA PARTICION DEL ESPACIO (phi,u), COMENZANDO A LLENAR 
C      COLUMNA POR COLUMNA.
C      PARTICION DE PHI       
       do i=1,int(fa1)
C      PARTICION DE  U      
         do j=1,int(fa2)
C      A CADA BIN (i,j) SE LE ASIGNAN 0 ELEMENTOS DENTRO DE SI          
            his(i,j)=0.d0
         enddo
       enddo
C      EN ESTE PASO SE CALCULA EL NUMERO DE ELEMENTOS EN CADA BIN       
       do l=1,ndata
         i=int(x(l)*fa1)+1
         j=int(y(l)*fa2)+1
         his(i,j)=his(i,j)+1.
       enddo
C      CALCULO DE LA ENTROPIA DE INFORMACION
C      BARRIENDO COLUMNA POR COLUMNA       
       s=0.
       do i=1,int(fa1)
         do j=1,int(fa2)
            prob=his(i,j)/dfloat(ndata)
            if(his(i,j).gt.0) s = s - prob * log(prob)
         enddo
       enddo
C      ENTROPIA NORMALIZADA       
       s=s/dlog(fa1*fa2)
       return
       END

