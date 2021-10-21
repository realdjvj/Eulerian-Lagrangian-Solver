!        6  2  5
!         \ | /
!        3--0--1
!         / | \
!        7  4  8

PROGRAM MAIN
    IMPLICIT NONE
    INTEGER,PARAMETER   ::  dp = selected_real_kind(15,307)
    REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  U,V,rho,U2,V2,rho2
    REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE  ::  f,fstar,f_eq
    REAL(kind=dp),DIMENSION(:)  ::  w(0:8)
    INTEGER,DIMENSION(:)        ::  cx(0:8),cy(0:8)
    INTEGER             ::  N,Nx,Ny,i,j,k,t,dt,T_end
    REAL(kind=dp)       ::  U0,tau,omega,U_res1,V_res1,r_res1,U_res_norm,V_res_norm,r_res_norm
    REAL(kind=dp)       ::  jet_height
    !!!!! Open Inputs !!!!! 
    OPEN(1,FILE='input.txt',status='old')
    READ(1,*)
    READ(1,*)   N,U0,dt,tau,jet_height,T_end
    CLOSE(1)
    Nx = N 
    Ny = N 
    CALL init(U,V,rho,f,Nx,Ny,U0)
    ALLOCATE(f_eq(0:Nx,0:Ny,0:8),fstar(0:Nx,0:Ny,0:8))
    ALLOCATE(U2(0:Nx,0:Ny),V2(0:Nx,0:Ny),rho2(0:Nx,0:Ny))
    w(0) = 4.0_dp/9.0_dp 
    w(1:4) = 1.0_dp/9.0_dp 
    w(4:) = 1.0_dp/36.0_dp
    cx = (/0,1,0,-1,0,1,-1,-1,1/)
    cy = (/0,0,1,0,-1,1,1,-1,-1/)
    omega = dt/tau
    t = 0
    U_res_norm = 1.0
    V_res_norm = 1.0
    r_res_norm = 1.0
    DO WHILE(t<T_end)! &
        ! .AND. U_res_norm > 1e-6 &
        ! .AND. V_res_norm > 1e-6 & 
        ! .AND. r_res_norm > 1e-6)
        !!!!! Collision !!!!!
        CALL f_eq_calc(U,V,rho,f_eq) 
        fstar = (1.0_dp - omega)*f + omega*f_eq
        !!!!! Streaming !!!!! 
        ! DO i = 1,Nx-1
        !     DO j = 1,Ny-1
        !         f(i+1,j,1) = fstar(i,j,1)
        !         f(i,j-1,2) = fstar(i,j,2)
        !         f(i-1,j,3) = fstar(i,j,3)
        !         f(i,j+1,4) = fstar(i,j,4)
        !         f(i+1,j-1,5) = fstar(i,j,5)
        !         f(i-1,j-1,6) = fstar(i,j,6)
        !         f(i-1,j+1,7) = fstar(i,j,7)
        !         f(i+1,j+1,8) = fstar(i,j,8)
        !     END DO 
        ! END DO
        !!!!! Matrix indexed streaming !!!!!
        f(:,:,0)            = fstar(:,:,0)
        f(2:Nx,1:Ny-1,1)    = fstar(1:Nx-1,1:Ny-1,1)
        f(1:Nx-1,0:Ny-2,2)  = fstar(1:Nx-1,1:Ny-1,2)
        f(0:Nx-2,1:Ny-1,3)  = fstar(1:Nx-1,1:Ny-1,3)
        f(1:Nx-1,2:Ny,4)    = fstar(1:Nx-1,1:Ny-1,4)
        f(2:Nx,0:Ny-2,5)    = fstar(1:Nx-1,1:Ny-1,5)
        f(0:Nx-2,0:Ny-2,6)  = fstar(1:Nx-1,1:Ny-1,6)
        f(0:Nx-2,2:Ny,7)    = fstar(1:Nx-1,1:Ny-1,7)
        f(2:Ny,2:Ny,8)      = fstar(1:Nx-1,1:Ny-1,8)

        CALL bcs_jet(f,fstar,Nx,Ny,U0,jet_height,w)
        CALL reconstruct(f,U2,V2,rho2,cx,cy)
        IF(MOD(t,1000) == 0)THEN 
            PRINT *,'========================================='
            PRINT 101,t,N,U0
            101 FORMAT(' t =',I10,'   N =',I5,'   U0 =',F8.3)
            PRINT 201,tau
            201 FORMAT(' Tau =',F16.8)
            PRINT *,'- - - - - - - - - - - - - - - - - - - - -'
            PRINT 102,MAXVAL(U)
            102 FORMAT(' max U vel      =',EN17.8)
            PRINT 103,MAXVAL(V)
            103 FORMAT(' max V vel      =',EN17.8)
            PRINT 104,r_res_norm
            104 FORMAT(' res Rho        =',EN17.8)
            PRINT 105,U_res_norm
            105 FORMAT(' res U vel      =',EN17.8)
            PRINT 106,V_res_norm
            106 FORMAT(' res V vel      =',EN17.8)
        END IF
        t = t+dt
        IF(t == dt)THEN
            U_res1 = NORM2(U2-U)
            V_res1 = NORM2(V2-V)
            r_res1 = NORM2(rho2-rho)
        END IF
        U_res_norm = NORM2(U2-U)/U_res1                
        V_res_norm = NORM2(V2-V)/V_res1        
        r_res_norm = NORM2(rho2-rho)/r_res1
        U(:,:) = U2(:,:)
        V(:,:) = V2(:,:)
        rho(:,:) = rho2(:,:)
    END DO


    CALL output(U,V,rho,U0)
    DEALLOCATE(U,V,rho,f,f_eq,fstar)
    CONTAINS
        SUBROUTINE init(U,V,rho,f,Nx,Ny,U0)
            IMPLICIT NONE
            INTEGER     ::  Nx,Ny
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  U,V,rho
            REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE  ::  f
            REAL(kind=dp)   ::  U0 
            ALLOCATE(U(0:Nx,0:Ny),V(0:Nx,0:Ny),rho(0:Nx,0:Ny))
            ALLOCATE(f(0:Nx,0:Ny,0:8))
            U(:,:) = 0.0
            V(:,:) = 0.0
            rho(:,:) = 1.0_dp
            CALL f_eq_calc(U,V,rho,f)
            RETURN
        END SUBROUTINE init

        SUBROUTINE bcs_jet(f,fstar,Nx,Ny,U0,jet_height,w)
            IMPLICIT NONE
            INTEGER ::  Nx,Ny
            REAL(kind=dp)   ::  U0,jet_height
            INTEGER         ::  jet_top,jet_bot
            REAL(kind=dp),DIMENSION(:)  ::  w(0:8)
            REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE  ::  f,fstar
            REAL(kind=dp)      ::  rho_in
            ! bounce back in the jet inlet
            ! open boundaries everywhere else
            jet_top = int((Ny-jet_height)/2.0)
            jet_bot = int((Ny+jet_height)/2.0)
            rho_in  = SUM(f)/SIZE(f)*9.0_dp
            f(1,jet_top:jet_bot,1) = fstar(1,jet_top:jet_bot,3) - w(3)*cx(3)*rho_in*U0*3.0_dp
            f(1,jet_top:jet_bot,5) = fstar(1,jet_top:jet_bot,7) - w(7)*cx(7)*rho_in*U0*3.0_dp
            f(1,jet_top:jet_bot,8) = fstar(1,jet_top:jet_bot,6) - w(6)*cx(6)*rho_in*U0*3.0_dp

            !left walls
            f(1,0:jet_top,1) = fstar(1,0:jet_top,3)
            f(1,0:jet_top,5) = fstar(1,0:jet_top,7)
            f(1,0:jet_top,8) = fstar(1,0:jet_top,6)
            f(1,jet_bot:Ny,1) = fstar(1,jet_bot:Ny,3)
            f(1,jet_bot:Ny,5) = fstar(1,jet_bot:Ny,7)
            f(1,jet_bot:Ny,8) = fstar(1,jet_bot:Ny,6)

            !OPEN BOUNDARY
            !right exit
            f(Nx-1,:,0)         = fstar(Nx-1,:,0)
            f(Nx-1,:,1)         = fstar(Nx,:,1)
            f(Nx-1,1:Ny,2)      = fstar(Nx-1,0:Ny-1,2)
            f(Nx-1,:,3)         = fstar(Nx-2,:,3)
            f(Nx-1,0:Ny-1,4)    = fstar(Nx-1,1:Ny,4)
            f(Nx-1,1:Ny,5)      = fstar(Nx,0:Ny-1,5)
            f(Nx-1,1:Ny,6)      = fstar(Nx-2,0:Ny-1,6)
            f(Nx-1,0:Ny-1,7)    = fstar(Nx-2,1:Ny,7)
            f(Nx-1,0:Ny-1,8)    = fstar(Nx,1:Ny,8)

            !top exit
            f(:,1,0)        = fstar(:,1,0)
            f(0:Nx-1,1,1)   = fstar(1:Nx,1,1)
            f(:,1,2)        = fstar(:,0,2)
            f(1:Nx,1,3)     = fstar(0:Nx-1,1,3)
            f(:,1,4)        = fstar(:,2,4)
            f(0:Nx-1,1,5)   = fstar(1:Nx,0,5)
            f(1:Nx,1,6)     = fstar(0:Nx-1,0,6)
            f(1:Nx,1,7)     = fstar(0:Nx-1,2,7)
            f(0:Nx-1,1,8)   = fstar(1:Nx,2,8)

            !bottom exit
            f(:,Ny-1,0)        = fstar(:,Ny-1,0)
            f(0:Nx-1,Ny-1,1)   = fstar(1:Nx,Ny-1,1)
            f(:,Ny-1,2)        = fstar(:,Ny,2)
            f(1:Nx,Ny-1,3)     = fstar(0:Nx-1,Ny-1,3)
            f(:,Ny-1,4)        = fstar(:,Ny-2,4)
            f(0:Nx-1,Ny-1,5)   = fstar(1:Nx,Ny,5)
            f(1:Nx,Ny-1,6)     = fstar(0:Nx-1,Ny,6)
            f(1:Nx,Ny-1,7)     = fstar(0:Nx-1,Ny-2,7)
            f(0:Nx-1,Ny-1,8)   = fstar(1:Nx,Ny-2,8)

        END SUBROUTINE bcs_jet

        SUBROUTINE bcs_lid(f,fstar,Nx,Ny,U0)
            IMPLICIT NONE
            INTEGER ::  Nx,Ny
            REAL(kind=dp)   ::  U0
            REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE  ::  f,fstar 
            REAL(kind=dp),DIMENSION(:),ALLOCATABLE      ::  rho_top 
            ALLOCATE(rho_top(0:Nx))
            !!!!! Left wall !!!!! 
            f(1,1:Ny-2,1) = fstar(1,1:Ny-2,3)
            f(1,1:Ny-2,5) = fstar(1,1:Ny-2,7)
            f(1,1:Ny-2,8) = fstar(1,1:Ny-2,6)
            !!!!! Bottom wall !!!!!
            f(1:Nx-2,Ny,2) = fstar(1:Nx-2,Ny,4)
            f(1:Nx-2,Ny,6) = fstar(1:Nx-2,Ny,8)
            f(1:Nx-2,Ny,5) = fstar(1:Nx-2,Ny,7)
            !!!!! Right wall !!!!!
            f(Nx,1:Ny-2,3) = fstar(Nx,1:Ny-2,1)
            f(Nx,1:Ny-2,6) = fstar(Nx,1:Ny-2,8)
            f(Nx,1:Ny-2,7) = fstar(Nx,1:Ny-2,5)
            !!!!! Top wall (moving) !!!!! 
            rho_top = fstar(:,Ny,0) + fstar(:,Ny,1) + fstar(:,Ny,3) + 2.0_dp*(fstar(:,Ny,4) + fstar(:,Ny,7) + fstar(:,Ny,8))
            f(1:Nx-2,1,4) = fstar(1:Nx-2,1,2) 
            f(1:Nx-2,1,7) = fstar(1:Nx-2,1,5) - rho_top(1:Nx-2) * U0 / 6.0_dp
            f(1:Nx-2,1,8) = fstar(1:Nx-2,1,6) + rho_top(1:Nx-2) * U0 / 6.0_dp

            !!!!! CORNERS !!!!!
            !!!!! Top left !!!!!
            f(1,1,1) = fstar(1,1,3)
            f(1,1,8) = fstar(1,1,6)
            f(1,1,4) = fstar(1,1,2)
            !!!!! bottom left !!!!!
            f(1,Ny-1,1) = fstar(1,Ny-1,3)
            f(1,Ny-1,5) = fstar(1,Ny-1,7)
            f(1,Ny-1,2) = fstar(1,Ny-1,4)
            !!!!! top right !!!!!               
            f(Nx-1,1,3) = fstar(Nx-1,1,1)               
            f(Nx-1,1,7) = fstar(Nx-1,1,5)               
            f(Nx-1,1,4) = fstar(Nx-1,1,2)
            !!!!! bottom right !!!!!
            f(Nx-1,Ny-1,3) = fstar(Nx-1,Ny-1,1)
            f(Nx-1,Ny-1,6) = fstar(Nx-1,Ny-1,8)
            f(Nx-1,Ny-1,2) = fstar(Nx-1,Ny-1,4)
            DEALLOCATE(rho_top)
            RETURN 
        END SUBROUTINE bcs_lid

        SUBROUTINE f_eq_calc(U,V,rho,f_eq)
            IMPLICIT NONE
            ! INTEGER ::  i 
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE   ::  U,V,rho
            REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE  ::  f_eq
                f_eq(:,:,0) = 2.0_dp*rho/9.0_dp * (2.0_dp - 3.0_dp * (U**2 + V**2))
                f_eq(:,:,1) = rho/18.0_dp * (2.0_dp + 6.0_dp*U + 9.0_dp*U**2 - 3.0_dp * (U**2 + V**2))
                f_eq(:,:,2) = rho/18.0_dp * (2.0_dp + 6.0_dp*V + 9.0_dp*V**2 - 3.0_dp * (U**2 + V**2))
                f_eq(:,:,3) = rho/18.0_dp * (2.0_dp - 6.0_dp*U + 9.0_dp*U**2 - 3.0_dp * (U**2 + V**2))
                f_eq(:,:,4) = rho/18.0_dp * (2.0_dp - 6.0_dp*V + 9.0_dp*V**2 - 3.0_dp * (U**2 + V**2))
                f_eq(:,:,5) = rho/36.0_dp * (1.0_dp + 3.0_dp*(U+V) + 9.0_dp*U*V + 3.0_dp*(U**2 + V**2))
                f_eq(:,:,6) = rho/36.0_dp * (1.0_dp - 3.0_dp*(U-V) - 9.0_dp*U*V + 3.0_dp*(U**2 + V**2))
                f_eq(:,:,7) = rho/36.0_dp * (1.0_dp - 3.0_dp*(U+V) + 9.0_dp*U*V + 3.0_dp*(U**2 + V**2))
                f_eq(:,:,8) = rho/36.0_dp * (1.0_dp + 3.0_dp*(U-V) - 9.0_dp*U*V + 3.0_dp*(U**2 + V**2))
            RETURN 
        END SUBROUTINE f_eq_calc

        SUBROUTINE reconstruct(f,U,V,rho,cx,cy)
            IMPLICIT NONE
            INTEGER         ::  i,j
            INTEGER,DIMENSION(:)        ::  cx(0:8),cy(0:8)        
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  U,V,rho 
            REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE  ::  f 
            DO i = 0,ubound(U,dim=1)
                DO j = 0,ubound(U,dim=1)
                    rho(i,j) = SUM(f(i,j,:))
                    U(i,j) = ((f(i,j,1)+f(i,j,5)+f(i,j,8)) - (f(i,j,3)+f(i,j,6)+f(i,j,7)))/rho(i,j)
                    V(i,j) = ((f(i,j,2)+f(i,j,5)+f(i,j,6)) - (f(i,j,4)+f(i,j,7)+f(i,j,8)))/rho(i,j)
                END DO 
            END DO 
            RETURN 
        END SUBROUTINE reconstruct

        SUBROUTINE output(U,V,rho,U0,t)
            IMPLICIT NONE
            INTEGER ::  i,j 
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  U,V,rho 
            REAL(kind=dp)   ::  U0,t
            OPEN(3,FILE='U.csv',STATUS='replace')
            OPEN(4,FILE='V.csv',STATUS='replace')
            OPEN(5,FILE='rho.csv',STATUS='replace')
            DO j = 1,ubound(U,dim=2)-1
                WRITE(3,101) (U(i,j)/U0, i = 1,ubound(U,dim=1)-1)
                WRITE(4,101) (V(i,j)/U0, i = 1,ubound(V,dim=1)-1)
                WRITE(5,101) (rho(i,j), i = 1,ubound(rho,dim=1)-1)
                101 FORMAT(*(E16.8))
            END DO 
            CLOSE(3)
            CLOSE(4)
            CLOSE(5)
            RETURN
        END SUBROUTINE output
END PROGRAM MAIN