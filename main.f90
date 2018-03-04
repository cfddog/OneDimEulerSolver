!//modules for contain the variables    
    module MD_VARIABLES
      real,pointer,dimension(:,:) :: Qnew,Qold,QTMP,DQ
      real,pointer,dimension(:)   :: den,vel,pre
      real,pointer,dimension(:,:) :: fp,fn,flux,HFLUX
      REAL,POINTER,DIMENSION(:)   :: X 
      !//DIAG MATRIX
      REAL,POINTER,DIMENSION(:,:) :: ARR
      REAL,POINTER,DIMENSION(:  ) :: BRR
      REAL,POINTER,DIMENSION(:  ) :: XRR
    end module
!//module for containing the control parameters
    module MD_CTRL
      integer Np
      integer vL
      integer itype_ischm,itype_split
      real    dt,Nt
      REAL :: GAMA=1.0
      REAL :: PAI=3.14159265
    end module
!//main program     
    program main
      !*********************************************
      !INITIALIZATION
      !*********************************************
      !the number of points and the length of domain
      Np=400
      MP=3
      vL=10.
      DT=0.001
      NT=1800
      !allocation
      ALLOCATE(X(Np))
      allocate(Qnew(Np,3))
      allocate(Qold(Np,3))
      ALLOCATE(QTMP(Np,3),DQ(Np,3))
      allocate(den(Np),vel(Np),pre(Np))
      allocate(fp(Np,3),fn(Np,3),flux(Np,3))
      ALLOCATE(HFLUX(0:Np,3)) !RECONSTRUCT FLUX CONTAINING TWO EXTRAPOLATION POINTS
      ALLOCATE(ARR(0:Np,MP),BRR(0:NP),XRR(0:NP))
      !INIT THE FIELD
      DX=VL/(NP-1)
      DO I=1,400
         X(I)=(I-1)*DX
         IF(X(I) .LT. 1.) THEN
             DEN(I)=3.857143
             VEL(I)=2.629369
             PRE(I)=10.33333
         ELSE
             DEN(I)=1+0.2*SIN(5.*X(I))
             VEL(I)=0.0
             PRE(I)=1.0
         ENDIF
      ENDDO
      DO I=1,400
          QNEW(I,1)=DEN(I)
          QNEW(I,2)=DEN(I)*VEL(I)
          QNEW(I,3)=DEN(I)*(0.5*VEL(I)**2+PRE(I)/DEN(I)/(GAMA-1))
          QOLD(I,:)=QNEW(I,:)
      ENDDO
      !TIME ADVANCE
       DO IT=1,NT
          CALL TIME_ADVANCE()
       ENDDO   
      ! 
    stop
    end program
    
    SUBROUTINE TIME_ADVANCE()
    USE md_variables
    USE md_ctrl
    !FOUR STEP RUNGE-KUTTA METHOD
    !STEP ONE
     CALL RHS()
     DO IP=1,NP
        DO IV=1,3 
           QNEW(IP,IV)=QOLD(IP,IV)+DQ(IP,IV)*DT/2.
           QTMP(IP,IV)=QNEW(IP,IV)
        ENDDO
        DEN(IP)=QNEW(IP,1)
        VEL(IP)=QNEW(IP,2)/QNEW(IP,1)
        PRE(IP)=(GAMA-1)*(QNEW(IP,3)-DEN(I)*0.5*VEL(I)**2)
     ENDDO
    !STEP TWO
     CALL RHS()
     DO IP=1,NP
        DO IV=1,3 
           QNEW(IP,IV)=QOLD(IP,IV)+DQ(IP,IV)*DT/2.
           QTMP(IP,IV)=QTMP(IP,IV)+2.*QNEW(IP,IV)
        ENDDO
        DEN(IP)=QNEW(IP,1)
        VEL(IP)=QNEW(IP,2)/QNEW(IP,1)
        PRE(IP)=(GAMA-1)*(QNEW(IP,3)-DEN(I)*0.5*VEL(I)**2)
     ENDDO
     !STEP THREE
     CALL RHS()
     DO IP=1,NP
        DO IV=1,3 
           QNEW(IP,IV)=QOLD(IP,IV)+DQ(IP,IV)*DT
           QTMP(IP,IV)=QTMP(IP,IV)+QNEW(IP,IV)
        ENDDO
        DEN(IP)=QNEW(IP,1)
        VEL(IP)=QNEW(IP,2)/QNEW(IP,1)
        PRE(IP)=(GAMA-1)*(QNEW(IP,3)-DEN(I)*0.5*VEL(I)**2)
     ENDDO 
     !STEP FOUR
     CALL RHS()
     DO IP=1,NP
        DO IV=1,3 
           QNEW(IP,IV)=(QTMP(IP,IV)-QOLD(IP,IV))/3.+DQ(IP,IV)*DT/6.
        ENDDO
        DEN(IP)=QNEW(IP,1)
        VEL(IP)=QNEW(IP,2)/QNEW(IP,1)
        PRE(IP)=(GAMA-1)*(QNEW(IP,3)-DEN(I)*0.5*VEL(I)**2)
     ENDDO
     QOLD=QNEW
     QTMP=0.0
     
    RETURN
    END SUBROUTINE
    
    SUBROUTINE RHS()
    USE md_variables
    USE md_ctrl
    !SPILT METHOD
    CALL FLUXSPLIT
    CALL FLUXRECONSTRUCT
    RETURN
    END SUBROUTINE
    
    SUBROUTINE FLUXSPLIT()
    USE md_variables
    USE md_ctrl
    REAL VLAMDA(3),VLAMDA_P(3),VLAMDA_N(3)
    REAL :: EPS=1.E-6
    !//FLUX AT THE POINTS
    DO IP=1,NP
       FLUX(IP,1)=DEN(IP)
       FLUX(IP,2)=DEN(IP)*VEL(IP)**2+PRE(IP)
       FLUX(IP,3)=VEL(IP)*(PRE(IP)+QNEW(IP,3))
       !//POSITIVE SPLIT
       AVEL=SQRT(GAMA*PRE(IP)/DEN(IP))
       VLAMDA(1)=VEL(IP)
       VLAMDA(2)=VEL(IP)-AVEL
       VLAMDA(3)=VEL(IP)+AVEL
       DO I1=1,3
          VLAMDA_P(I1)=(VLAMDA(I1)+SQRT(VLAMDA(I1)**2+EPS**2))/2.
       ENDDO
       W=(3-GAMA)*(VLAMDA_P(2)+VLAMDA_P(3))*AVEL**2/2./(GAMA-1)
       FP(IP,1)=( 2*(GAMA-1)*VLAMDA_P(1)+VLAMDA_P(2)+VLAMDA_P(3) )*DEN(IP)/2./GAMA
       FP(IP,2)=( 2*(GAMA-1)*VLAMDA_P(1)*VEL(IP)+VLAMDA_P(2)*(VEL(IP)-AVEL)+VLAMDA_P(3)*(VEL(IP)+AVEL))*DEN(IP)/2./GAMA
       FP(IP,3)=(   (GAMA-1)*VLAMDA_P(1)*VEL(IP)**2+VLAMDA_P(2)*(VEL(IP)-AVEL)**2/2.+VLAMDA_P(3)*(VEL(IP)+AVEL)**2/2.+W)*DEN(IP)/2./GAMA
       !NEGATIVE SPLIT
       DO I1=1,3
          VLAMDA_N(I1)=(VLAMDA(I1)-SQRT(VLAMDA(I1)**2+EPS**2))/2.
       ENDDO
       W=(3-GAMA)*(VLAMDA_N(2)+VLAMDA_N(3))*AVEL**2/2./(GAMA-1)
       FN(IP,1)=( 2*(GAMA-1)*VLAMDA_N(1)+VLAMDA_N(2)+VLAMDA_N(3) )*DEN(IP)/2./GAMA
       FN(IP,2)=( 2*(GAMA-1)*VLAMDA_N(1)*VEL(IP)+VLAMDA_N(2)*(VEL(IP)-AVEL)+VLAMDA_N(3)*(VEL(IP)+AVEL))*DEN(IP)/2./GAMA
       FN(IP,3)=(   (GAMA-1)*VLAMDA_N(1)*VEL(IP)**2+VLAMDA_N(2)*(VEL(IP)-AVEL)**2/2.+VLAMDA_N(3)*(VEL(IP)+AVEL)**2/2.+W)*DEN(IP)/2./GAMA
    ENDDO
    
    RETURN
    END SUBROUTINE
    
    SUBROUTINE FLUXRECONSTRUCT
    USE MD_VARIABLES
    USE MD_CTRL
    !//POSITIVE
    DO IP=0,NP
        !//SPECIFY THE ARRAY
        ARR=0.0
        BRR=0.0
        XRR=0.0
        IF(IP .EQ. 0) THEN
        ARR(IP,1)=0;ARR(IP,2)=1.;ARR(IP,3)=0
        BRR(IP)
        ELSEIF(IP .EQ. 1) THEN
        ARR(IP,1)=0;ARR(IP,2)=1,ARR(IP,3)=0.
        ELSEIF(IP .EQ. NP) THEN
        ARR(IP,2)=1.    
        ELSE
        ARR(IP,1)=3.;ARR(IP,2)=6.;ARR(IP,3)=1.    
        ENDIF
        
    ENDDO
    
    RETURN
    END SUBROUTINE