PROGRAM main
    USE fhdict
    IMPLICIT NONE
    INTEGER,PARAMETER   ::  dp = selected_real_kind(15,307)
    REAL(kind=dp)   ::  Re, dx, dy, H, dt, P_inter,P_res1,P_res_norm,U_res_norm,V_res_norm
    REAL(kind=dp)   ::  r1,r2,r3,r4,Q,U_res1,V_res1,tol,P_res,U_res,V_res

    REAL(kind=dp)   ::  Txxe,Txxw,Txyn,Txys,uee,uww,unn,uss,vnn,vss
    REAL(kind=dp)   ::  Tyyn,Tyys,Txye,Txyw,vee,vww

    REAL(kind=dp)   ::  res_normalized,res1,t,U0,Tend,L
    REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  P,U,V,U2,V2,P2,res,velmag,Eg,Ep,Fdx,Fdy
    INTEGER         ::  Nx,Ny,i,j,k,jet_height,N,t_count

    TYPE(hashtable)  ::  particles
    !!!!! READ INPUTS !!!!!
    OPEN(1,FILE='input.txt',STATUS='old')
    READ(1,*)
    READ(1,*) Re, N, L, U0, dt, jet_height, Tend
    CLOSE(1)
    Nx = N
    Ny = N
    ALLOCATE(Eg(0:Nx+1,0:Ny+1),Ep(0:Nx+1,0:Ny+1))
    ALLOCATE(Fdx(0:Nx+1,0:Ny+1),Fdy(0:Nx+1,0:Ny+1))
    Eg(:,:) = 1.0_dp
    Ep(:,:) = 0.0_dp
    Fdx(:,:) = 0.0_dp
    Fdy(:,:) = 0.0_dp
    !!!!! INITIALIZE !!!!!!
    CALL particles%init(1000000)
    CALL particles%put(key=1,rvals=[0.0001_dp, 0.0_dp, 0.5_dp, U0, 0.0_dp])
    CALL seed_particles(particles,U0)
    dx = L/Nx
    dy = L/Ny
    CALL init(P,U,V,P2,U2,V2,Nx,Ny,U0)
    CALL bcs(P,U,V,jet_height,U0)
    ALLOCATE(res(0:Nx+1,0:Ny+1),velmag(0:Nx+1,0:Ny+1))
    t = 0.0_dp
    t_count = 0
    P_res_norm = 1.0_dp
    U_res_norm = 1.0_dp
    V_res_norm = 1.0_dp
    r1 = 1.0_dp/(dx)
    r2 = 1.0_dp/(dy)
    r3 = dt/(dx**2)
    r4 = dt/(dy**2)
    !!!!! TIME LOOP !!!!!
    DO WHILE(t < Tend)
        !!!!! CONVECTIVE + VISCOUS TERMS !!!!!
        DO i = 1,ubound(U,dim=1)-1
            DO j = 1,ubound(U,dim=2)-1
                Txxe = 2.0_dp/Re*(U(i+1,j) - U(i,j))/dx
                Txxw = 2.0_dp/Re*(U(i,j) - U(i-1,j))/dx
                Txyn = ((V(i+1,j-1) - V(i,j-1))/dx + (U(i,j-1)-U(i,j))/dy)/Re
                Txys = ((V(i+1,j) - V(i,j))/dx + (U(i,j)-U(i,j+1))/dy)/Re
                uee = 0.5_dp*(U(i,j)+U(i+1,j))
                uww = 0.5_dp*(U(i,j)+U(i-1,j))
                unn = 0.5_dp*(U(i,j)+U(i,j-1))
                uss = 0.5_dp*(U(i,j)+U(i,j+1))
                vnn = 0.5_dp*(V(i,j-1)+V(i+1,j-1))
                vss = 0.5_dp*(V(i,j)+V(i+1,j))
                H = ((Txxe-Txxw)/dx + (Txyn-Txys)/dy) - (uee**2-uww**2)/dx - (unn*vnn-uss*vss)/dy &
                    - ((Fdx(i,j) + Fdx(i+1,j))/2.0_dp)
                U2(i,j) = U(i,j) + dt*H

                Tyyn = 2.0_dp/Re*(V(j,i-1)-V(j,i))/dy
                Tyys = 2.0_dp/Re*(V(j,i)-V(j,i+1))/dy
                Txye = ((V(j+1,i)-V(j,i))/dx + (U(j,i) - U(j,i+1))/dy)/Re
                Txyw = ((V(j,i)-V(j-1,i))/dx + (U(j-1,i)-U(j-1,i+1))/dy)/Re
                vnn = 0.5_dp*(V(j,i)+V(j,i-1))
                vss = 0.5_dp*(V(j,i)+V(j,i+1))
                uee = 0.5_dp*(U(j,i)+U(j,i+1))
                uww = 0.5_dp*(U(j-1,i)+U(j-1,i+1))
                vee = 0.5_dp*(V(j,i)+V(j+1,i))
                vww = 0.5_dp*(V(j,i)+V(j-1,i))
                H = (Tyyn-Tyys)/dy + (Txye-Txyw)/dx - (vnn**2 - vss**2)/dy - (uee*vee - uww*vww)/dx &
                    - ((Fdy(i,j) + Fdy(i,j+1))/2.0_dp)
                V2(j,i) = V(j,i) + dt*H
            END DO
        END DO

        CALL bcs(P,U2,V2,jet_height,U0)

        !!!!! POISSON'S EQUATIONS !!!!!
        k = 0
        res_normalized = 1
        P2(:,:) = P(:,:)
        DO WHILE(res_normalized > 1e-6_dp .AND. k < 10000)
            DO i = 1,ubound(P,dim=1)-1
                DO j = 1,ubound(P,dim=2)-1
                    Q = (U2(i,j)/dx - U2(i-1,j)/dx + V2(i,j-1)/dy - V2(i,j)/dy)/dt
                    P_inter = (Q - P(i+1,j)/dx**2 - P(i-1,j)/dx**2 - P(i,j-1)/dy**2 - P(i,j+1)/dy**2) &
                        /(-2.0_dp/dx**2 - 2.0_dp/dy**2)
                    res(i,j) = P_inter - P(i,j)
                    P(i,j) = P_inter
                END DO
            END DO
            IF(k==0) THEN
                res1 = NORM2(res)
            END IF
            res_normalized = NORM2(res)/res1
            k = k+1
        END DO

        CALL bcs(P,U2,V2,jet_height,U0)

        IF(t_count==0) THEN
            U_res1 = NORM2(U2-U)
            P_res1 = NORM2(P2-P)
            V_res1 = NORM2(V2-V)
        END IF
        U_res = NORM2(U2-U)
        P_res = NORM2(P2-P)
        V_res = NORM2(V2-V)
        U_res_norm = U_res/U_res1
        P_res_norm = P_res/P_res1
        V_res_norm = V_res/V_res1
        DO i=1,ubound(U,dim=1)-1
            DO j=1,ubound(U,dim=2)-1
                U(i,j) = U2(i,j) - dt*((P(i+1,j)-P(i,j))/(dx))
                V(j,i) = V2(j,i) - dt*((P(j,i)-P(j,i+1))/(dy))
            END DO
        END DO

        CALL bcs(P,U,V,jet_height,U0)
        CALL seed_particles(particles,U0)
        CALL move_particles(particles,P,U,V,Ep,Eg,Fdx,Fdy,dt,dx,dy)

        t = t_count*dt
        IF(MOD(t_count,100)==0) THEN
            PRINT *,'============================================='
            PRINT 100,t,Re,dx
            100 FORMAT(' t = ',ES14.6,'   Re = ',F8.1,'   dx = ',F6.4)
            PRINT 110,k
            110 FORMAT(' k = ',I5)
            PRINT *,'- - - - - - - - - - - - - - - - - - - - - - -'
            PRINT 101,P_res_norm
            101 FORMAT(' P_res_norm      = ',EN14.5)
            PRINT 102,U_res_norm
            102 FORMAT(' U_res_norm      = ',EN14.5)
            PRINT 103,V_res_norm
            103 FORMAT(' V_res_norm      = ',EN14.5)
            PRINT 104,particles%count
            104 FORMAT(' N particles     = ',I10)
            CALL collocate(P,U,V,Nx,Ny,particles,t_count)
        END IF
        t_count = t_count+1
    END DO

    CALL collocate(P,U,V,Nx,Ny,particles,t_count)
    DEALLOCATE(P,U,V,P2,U2,V2)
    DEALLOCATE(res,velmag)
    CONTAINS

        SUBROUTINE init(P,U,V,P2,U2,V2,Nx,Ny,U0)
            IMPLICIT NONE
            INTEGER     :: Nx,Ny
            REAL(kind=dp)   ::  U0
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  P,U,V,P2,U2,V2

            ALLOCATE(U(0:Nx,0:Ny+1),V(0:Nx+1,0:Ny),P(0:Nx+1,0:Ny+1),U2(0:Nx,0:Ny+1),V2(0:Nx+1,0:Ny),P2(0:Nx+1,0:Ny+1))
            P(:,:)  = 0.0_dp
            U(:,:)  = 0.0_dp
            V(:,:)  = 0.0_dp
            P2(:,:) = 0.0_dp
            U2(:,:) = 0.0_dp
            V2(:,:) = 0.0_dp
            RETURN
        END SUBROUTINE init

        SUBROUTINE seed_particles(particles,U0)
            IMPLICIT NONE
            TYPE(hashtable) ::  particles
            REAL(kind=dp)   ::  U0, randn
            INTEGER,ALLOCATABLE ::  mykeys(:)
            INTEGER         ::  i
            CALL particles%keys(mykeys)
            DO i = 1,10
                CALL random_number(randn)
                CALL particles%put(key=mykeys(ubound(mykeys,dim=1))+i, &
                        rvals=[0.0001_dp, 0.0_dp,0.45_dp+(randn*0.1_dp),U0,0.0_dp])
            END DO
            RETURN
        END SUBROUTINE

        SUBROUTINE bcs(P,U,V,jet_height,U0)
            IMPLICIT NONE
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  P,U,V
            REAL(kind=dp)   ::  U0
            INTEGER         ::  jet_height,jet_bot,jet_top
            INTEGER,DIMENSION(:),ALLOCATABLE    ::  mykeys

            jet_top = (ubound(U,dim=2) - jet_height)/2
            jet_bot = (ubound(U,dim=2) + jet_height)/2
            U(:,ubound(U,dim=2)) = U(:,ubound(U,dim=2)-1)
            U(:,0) = U(:,1)
            U(ubound(U,dim=1),:) = U(ubound(U,dim=1)-1,:)
            U(0,0:jet_top) = 0.0_dp!U(1,0:jet_top)
            U(0,jet_top:jet_bot) = U0
            U(0,jet_bot:ubound(U,dim=2)) = 0.0_dp!U(1,jet_bot:ubound(U,dim=2))

            V(0,:) = 0.0_dp!V(1,:)
            V(ubound(V,dim=1),:) = V(ubound(V,dim=1)-1,:)
            V(:,0) = V(:,1)
            ! V(:,ubound(V,dim=2)) = U0 ! Jet in crossflow
            V(:,ubound(V,dim=2)) = V(:,ubound(V,dim=2)-1) ! Laminar Jet

            P(:,0) = P(:,1)
            P(0,:) = P(1,:)
            P(ubound(P,dim=1),:) = P(ubound(P,dim=1)-1,:)
            P(:,ubound(P,dim=2)) = P(:,ubound(P,dim=2)-1)
            RETURN
        END SUBROUTINE bcs

        SUBROUTINE move_particles(particles,P,U,V,Ep,Eg,Fdx,Fdy,dt,dx,dy)
            IMPLICIT NONE
            TYPE(hashtable)     ::  particles
            INTEGER,DIMENSION(:),ALLOCATABLE    ::  mykeys
            REAL(kind=dp)       ::  u_new,v_new,x_new,y_new,dt,U0,dx,dy
            INTEGER     ::  i,this_key,ip,jp
            REAL(kind=dp),DIMENSION(:),ALLOCATABLE  ::  array_old
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  P ,Eg,Ep,U,V,Fdx,Fdy
            REAL(kind=dp)   ::  dpdx,dpdy,r,pi=4.0_dp*DATAN(1.0_dp)!,m=1.0_dp

            CALL particles%keys(mykeys)
            Ep(:,:) = 0.0_dp
            DO i = 1,ubound(mykeys,dim=1)
                this_key = mykeys(i)
                CALL particles%get(key=this_key, rvals=array_old)
                r = array_old(1)
                ip = NINT((array_old(2) + dx/2.0_dp)/dx)
                jp = NINT((1.0_dp-array_old(3) + dy/2.0_dp)/dy)
                Ep(ip,jp) = Ep(ip,jp) + pi*r**2/(dx*dy)
            END DO
            DO i = 1,ubound(mykeys,dim=1)
                this_key = mykeys(i)
                CALL particles%get(key=this_key,rvals=array_old)
                ip = NINT((array_old(2)+dx/2.0_dp)/dx)
                jp = NINT((1.0_dp-array_old(3) + dy/2.0_dp)/dy)
                dpdx = (P(ip+1,jp)-P(ip-1,jp))/(2.0_dp*dx)
                dpdy = (P(ip,jp+1)-P(ip,jp-1))/(2.0_dp*dy)
                x_new = array_old(2) + array_old(4)*dt
                y_new = array_old(3) + array_old(5)*dt

                Fdx(ip,jp) = fd(r,U(ip,jp),array_old(4),Ep(ip,jp))
                Fdy(ip,jp) = fd(r,V(ip,jp),array_old(5),Ep(ip,jp))
                !!!! Density of ALuminum used !!!!
                u_new = array_old(4) + dt*(1.0_dp/2700.0_dp * dpdx + Fdx(ip,jp))
                v_new = array_old(5) + dt*(1.0_dp/2700.0_dp * dpdy + Fdy(ip,jp))

                IF(array_old(2)>1.0_dp .OR. array_old(3)>1.0_dp &
                .OR. array_old(2)<0.0_dp .OR. array_old(3)<0.0_dp)THEN
                    CALL particles%del(this_key)
                ELSE
                    CALL particles%put(key=this_key, rvals=[array_old(1),x_new,y_new,u_new,v_new])
                END IF
            END DO
            Eg(:,:) = 1.0_dp-Ep(:,:)
            RETURN
        END SUBROUTINE move_particles

        SUBROUTINE drag_force(r,Ug,Up,Ep,fdx)
            IMPLICIT NONE
            REAL(kind=dp)   ::  fdx,r,Up,Ug,Ep,Eg,mu,rho,Vp,Re,Cd,beta,pi
            IF(Ep==0.0_dp)THEN
                fdx = 0.0_dp
            ELSE
                pi = 4.0_dp*DATAN(1.0_dp)
                rho = 1.225_dp
                mu = 1.5e-5_dp
                ! Eg = 1.0_dp-Ep
                Vp = pi*r**2
                Re = 2.0_dp*r*rho*ABS(Ug-Up)/mu
                IF(Re > 1000.0_dp)THEN
                    Cd = 0.44_dp
                ELSE
                    Cd = 24.0_dp/Re * (1.0_dp + 0.15*Re**(0.68_dp))
                ENDIF
                beta = 3.0_dp*Cd*Eg*Ep*rho*ABS(Ug-Up)/(4.0_dp*2.0_dp*r)
                fdx = beta*Vp/(Ep)*(Ug-Up)
            END IF
            print *,fdx
            RETURN
        END SUBROUTINE drag_force

        FUNCTION fd(r,Ug,Up,Ep_loc)
            IMPLICIT NONE
            REAL(kind=dp)   ::  fd,Vp,r,Re,rho,Ug,Up,mu,Ep_loc,beta,Cd,pii=4.0_dp*DATAN(1.0_dp)
            ! IF(Ep == 0.0_dp)THEN
            !     fd = 0.0_dp
            ! ELSE
                rho = 1.225_dp
                mu = 1.5e-5_dp
                Vp = pii*r**2
                Re = 2.0_dp*r*rho*ABS(Ug-Up)/mu
                IF(Re > 1000.0_dp)THEN
                    Cd = 0.44_dp
                ELSE
                    Cd = 24.0_dp/Re * (1.0_dp + 0.15_dp*Re**(0.68_dp))
                ENDIF
                ! Eg = 1.0_dp - Ep
                beta = 3.0_dp*Cd*(1.0_dp - Ep_loc)*Ep_loc*rho*ABS(Ug-Up)/(4.0_dp*2.0_dp*r)
                ! fd = beta*Vp/(Ep)*(Ug-Up)
                fd = beta/Ep_loc * (Ug - Up) / 2700.0_dp
            ! ENDIF
            RETURN
        END FUNCTION

        SUBROUTINE collocate(P,U,V,Nx,Ny,particles,t_count)
            IMPLICIT NONE
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  P,U,V
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  Pco,Uco,Vco
            INTEGER     ::  i,j,Nx,Ny,t_count
            TYPE(hashtable) ::  particles
            ALLOCATE(Pco(0:Nx,0:Ny),Uco(0:Nx,0:Ny),Vco(0:Nx,0:Ny))
            DO i = 0,Nx
                DO j = 0,Ny
                    Uco(i,j) = (U(i,j) + U(i,j+1))/2.0_dp
                    Vco(i,j) = (V(i,j) + V(i+1,j))/2.0_dp
                END DO
            END DO
            CALL output(P,Uco,Vco,particles,t_count)
            DEALLOCATE(Pco,Uco,Vco)
        END SUBROUTINE collocate

        SUBROUTINE output(P,U,V,particles,t_count)
            IMPLICIT NONE
            INTEGER     ::  i,j,t_count
            REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    ::  P,U,V
            INTEGER,DIMENSION(:),ALLOCATABLE            ::  mykeys
            INTEGER ::  this_key,count
            REAL(kind=dp),DIMENSION(:),ALLOCATABLE      ::  array
            TYPE(hashtable) ::  particles
            CHARACTER(len=10)   ::  end
            WRITE(end,"(I10)") t_count
            end = adjustl(end)
            CALL particles%keys(mykeys)
            OPEN(5,FILE='./data/particles'//trim(end)//'.csv',STATUS='replace')
            count = ubound(mykeys,dim=1)
            DO j = 1,count
                this_key = mykeys(j)
                CALL particles%get(key=this_key,rvals=array)
                WRITE(5,1) (array(i), i = 1,5)
            END DO
            CLOSE(5)

            OPEN(2,FILE='./data/U'//trim(end)//'.csv',STATUS='replace')
            DO j = 0,ubound(U,dim=2)
                WRITE(2,1) (U(i,j), i = 0,ubound(U,dim=1))
                1   FORMAT(*(F16.8))
           END DO
            CLOSE(2)
            OPEN(3,FILE='./data/P'//trim(end)//'.csv',STATUS='replace')
            DO j = 0,ubound(P,dim=2)
                WRITE(3,1) (P(i,j), i = 0,ubound(P,dim=1))
            END DO
            CLOSE(3)
            OPEN(4,FILE='./data/V'//trim(end)//'.csv',STATUS='replace')
            DO j = 0,ubound(V,dim=2)
                WRITE(4,1) (V(i,j), i = 0,ubound(V,dim=1))
            END DO
            CLOSE(4)
        END SUBROUTINE output

END PROGRAM main
