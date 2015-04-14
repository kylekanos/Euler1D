! Solves Eulerian Hydro equations in 1D

program Euler1D
    implicit none
    integer, parameter :: wp = selected_real_kind(15,307)
    integer, parameter :: nx = 100
    real(wp), parameter :: half = 0.5_wp
    real(wp), parameter :: zero = 0.0_wp
    
    real(wp), dimension(-2:nx+3):: ro, px, en, x
    real(wp) :: time, dt, dx, tend, CFL, gamma, gamma1
    integer :: step, i, nout

    dx=1.0/real(nx,wp)
    CFL = 0.3_wp
    gamma = 1.4
    gamma1 = gamma - 1.0_wp


! set up grid
    x(1) = half*dx
    do i=-2,nx+3
       x(i) = x(i-1) + dx
    enddo

    do i=-2,nx+3
       ro(i) = 1.16_wp
       px(i) = zero
       en(i) = 10.00_wp/gamma1
    enddo
   
    time = 0.0_wp
    tend = 0.2_wp
    step = 0
    nout = 0
    call output(x, ro, px, en, nout)
    
    print *," *** step *** | *** tnow *** | *** dt ***"
    do while(time < tend)
       
       dt=FindDT(ro, px, en)
       if (dt>(tend-time)) dt=(tend-time)
       step = step + 1
       print *,step, "|", time, "|", dt
       
       call AdvanceGrid(ro, px, en, x)
       time=time+dt
       
       if(mod(step,3) == 0) then
          nout = nout + 1
          call output(x, ro, px, en, nout)
       endif
    enddo
    print*, "endtime=",time

    call output(x, ro, px, en, nout+1)

 contains
   !> @brief computes new time-step
   real(wp) function FindDT(ro,px,e)
     real(wp), dimension(-2:nx+3), intent(in) :: ro,px,e
     real(wp) :: p,c,S,pom
     integer :: i
       
     S=0
     do i=1,nx
        p=(gamma-1._wp)*(en(i)-(1._wp/2._wp)*(px(i)*px(i))/ro(i))
        if (p<0) then
           write (*,*) "p<0!!!",en(i),ro(i),px(i)/ro(i)
           p=1e-10_wp
        endif
        c=SQRT(gamma*p/ro(i))
        pom=MAX(ABS(px(i)/ro(i)-c),ABS(px(i)/ro(i)+c))
        if(S<p) S=pom
     enddo
     FindDT=CFL*dx/S
   end function FindDT

   !> @brief writes the frame to curve file
   subroutine output(x,ro,px,e,n)
     real(wp),dimension(-2:):: x,ro,px,e
     integer i,n
     character(len=18) :: fname
     character(len=4) :: nchar
  
     write(nchar,'(i4.4)') n
     fname = "output"//nchar//".curve"
     open(unit=11,file=fname,status='replace')
     
     write(11,*) "# ro"
     do i=1,nx
        write(11,*) x(i), ro(i)
     enddo
     write(11,*) "# px"
     do i=1,nx
        write(11,*) x(i), px(i)
     enddo
     write(11,*) "# E"
     do i=1,nx
        write(11,*) x(i), en(i)
     enddo
     write(11,*) "# p"
     do i=1,nx
        write(11,*) x(i), (gamma-1._wp)*(en(i)-(1._wp/2._wp)*(px(i))*(px(i))/ro(i))
     enddo
     close(11)
      
     write (*,*) "Saved frame",n
   end subroutine output
   
   !> @brief marches the grid in time
   subroutine AdvanceGrid(ro,px,en,x)
     integer ::i
     real(wp), dimension(-2:nx+3) :: ro, px, en, x
     real(wp), dimension(-2:nx+3) :: ro2, px2, en2, dro, dpx, den
     
     call RKStep(dro,dpx,den,ro,px,en,x)
     ro2 = ro + dt*dro
     px2 = px + dt*dpx
     en2 = en + dt*den
     
     call RKStep(dro,dpx,den,ro2,px2,en2,x)
     ro = half*(ro + ro2 + dt*dro)
     px = half*(px + px2 + dt*dpx)
     en = half*(en + en2 + dt*den)
   end subroutine AdvanceGrid
   
   !> @brief RK step
   subroutine RKStep(dro,dpx,den,ro,px,e,x)
     real(wp), dimension(-2:nx+3) :: ro, px, e, dro, dpx, den, x
     real(wp), dimension(-2:nx+3) :: roL, pxL, enL, roR, pxR, enR
     
     dro=0
     dpx=0
     den=0
     
! set boundary condition
     call SetBC(ro,px,e) 
! set reconstruction
     call TVD_PWL(roL,pxL,enL,roR,pxR,enR,ro,px,e,x)
! determine fluxes
     !call CUMethod(dro,dpx,den,roL,pxL,enL,roR,pxR,enR)
     call KTMethod(dro,dpx,den,roL,pxL,enL,roR,pxR,enR,x)
   end subroutine RKStep
   
   !> @brief boundary conditions
   subroutine SetBC(ro,px,e)
     real(wp), dimension(-2:nx+3) :: ro, px, e  

     ro(-2:0) = ro(1)
     ro(nx+1:nx+3)=ro(nx)
      
     px(-2:0) = px(1)
     px(nx+1:nx+3)=-px(nx)

     en(-2:0)=10.001/(gamma-1._wp)
     en(nx+1:nx+3)=en(nx)
   end subroutine SetBC
   
   !> @brief Simple Piecewise-linear reconstruction
   subroutine TVD_PWL(roL,pxL,enL,roR,pxR,enR,ro,px,e,x)
     real(wp),dimension(-2:nx+3)::roL,pxL,enL,roR,pxR,enR,ro,px,e,x
     real(wp) :: sL,sC,sR
     integer :: i
     
     do i=0,nx+1
        sL = (ro(i  ) - ro(i-1))
        sR = (ro(i+1) - ro(i  ))
        sC = minmod(sL,sR)
        roR(i  ) = ro(i) - sC*half
        roL(i+1) = ro(i) + sC*half
     enddo
     do i=0,nx+1
        sL = (px(i  ) - px(i-1))       
        sR = (px(i+1) - px(i  ))
        sC = minmod(sL,sR)
        pxR(i  ) = px(i) - sC*half
        pxL(i+1) = px(i) + sC*half
     enddo
     do i=0,nx+1
        sL = (en(i  ) - en(i-1))
        sR = (en(i+1) - en(i  ))
        sC = minmod(sL,sR)
        enR(i  ) = en(i) - sC*half
        enL(i+1) = en(i) + sC*half
     enddo
   end subroutine TVD_PWL
   
   !> @brief MINMOD limiter
   real(wp) function minmod(a,b)
     real(wp), intent(in) :: a, b
     minmod = (sign(1._wp,a)+sign(1._wp,b))*min(abs(a),abs(b))/2._wp
   end function minmod
   
   !> @brief Computes Centered Upwind fluxes
   subroutine CUMethod(dro,dpx,den,roL,pxL,enL,roR,pxR,enR)
     real(wp), dimension(-2:nx+3) :: roL, pxL, enL, roR, pxR, enR
     real(wp), dimension(-2:nx+3) :: dro, dpx, den, aL, aR, pL, pR
     real(wp) :: cL, cR, p, dxi
     integer :: i
     
     dxi = 1.0_wp/dx
     
     do i=1,nx+1
        pL(i)=(gamma-1._wp)*(enL(i)-(1._wp/2._wp)*(pxL(i)*pxL(i))/roL(i))
        pR(i)=(gamma-1._wp)*(enR(i)-(1._wp/2._wp)*(pxR(i)*pxR(i))/roR(i))
        if ((pL(i) < zero).or.(pR(i) < zero)) then
           print*, "Error: negative pressure in CUMethod:",pL(i),pR(i)
           pL(i) = max(pL(i), 1e-10_wp)
           pR(i) = max(pR(i), 1e-10_wp)
        endif
        cL=SQRT(gamma*pL(i)/roL(i))
        cR=SQRT(gamma*pR(i)/roR(i))
        aL(i)=MIN(pxL(i)/roL(i) - cL, pxR(i)/roR(i) - cR, 0._wp)
        aR(i)=MAX(pxL(i)/roL(i) + cL, pxR(i)/roR(i) + cR, 0._wp)        
     enddo

! flux of ro: px
     dro(1)=-((aR(1)*pxL(1) - aL(1)*pxR(1))/(aR(1) - aL(1))&
             + (aL(1)*aR(1)/(aR(1) - aL(1)))*(roR(1) - roL(1)))*dxi
     do i=2,nx
        p = (aR(i)*pxL(i) - aL(i)*pxR(i))/(aR(i) - aL(i)) &
              + (aL(i)*aR(i)/(aR(i) - aL(i)))*(roR(i) - roL(i))
        dro(i-1) = dro(i-1) + p*dxi
        dro(i  ) = dro(i  ) - p*dxi
     enddo
     dro(nx)=dro(nx)+((aR(nx+1)*pxL(nx+1) - aL(nx+1)*pxR(nx+1))/(aR(nx+1) - aL(nx+1)) &
                       +(aL(nx+1)*aR(nx+1)/(aR(nx+1) - aL(nx+1)))*(roR(nx+1) - roL(nx+1)))*dxi
     dro=-dro

!flux of ro*u: ro*u*u + p
    dpx(1)=-((aR(1)*(pxL(1)**2/roL(1) + pL(1)) - aL(1)*(pxR(1)**2/roR(1) + pR(1)))/(aR(1)-aL(1))&
             + (aL(1)*aR(1)/(aR(1) - aL(1)))*(pxR(1) - pxL(1)))*dxi
     do i=2,nx
        p=(aR(i)*(pxL(i)**2/roL(i) + pL(i))&
               - aL(i)*(pxR(i)**2/roR(i) + pR(i)))/(aR(i) - aL(i))&
              +(aL(i)*aR(i)/(aR(i) - aL(i)))*(pxR(i) - pxL(i))
        dpx(i-1) = dpx(i-1) + p*dxi
        dpx(i  ) = dpx(i  ) - p*dxi
     enddo
     dpx(nx)=dpx(nx)+((aR(nx+1)*(pxL(nx+1)**2/roL(nx+1) + pL(nx+1))&
                       -aL(nx+1)*(pxR(nx+1)**2/roR(nx+1) + pR(nx+1)))/(aR(nx+1) - aL(nx+1))&
                       +(aL(nx+1)*aR(nx+1)/(aR(nx+1) - aL(nx+1)))*(pxR(nx+1) - pxL(nx+1)))*dxi
     dpx = -dpx
    
!flux of e: u*(e + p)
     den(1)=-((aR(1)*(pxL(1)/roL(1)*(enL(1)+pL(1)))&
             - aL(1)*(pxR(1)/roR(1)*(enR(1)+pR(1))))/(aR(1) - aL(1))&
             +(aL(1)*aR(1)/(aR(1) - aL(1)))*(enR(1) - enL(1)))*dxi
     do i=2,nx
        p=(aR(i)*(pxL(i)/roL(i)*(enL(i)+pL(i)))&
            -aL(i)*(pxR(i)/roR(i)*(enR(i)+pR(i))))/(aR(i)-aL(i))&
            +(aL(i)*aR(i)/(aR(i)-aL(i)))*(enR(i)-enL(i))
        den(i-1) = den(i-1) + p*dxi
        den(i  ) = den(i  ) - p*dxi
     enddo
     den(nx)=den(nx)+((aR(nx+1)*(pxL(nx+1)/roL(nx+1)*(enL(nx+1)+pL(nx+1)))&
                       -aL(nx+1)*(pxR(nx+1)/roR(nx+1)*(enR(nx+1)+pR(nx+1))))/(aR(nx+1)-aL(nx+1))&
                       +(aL(nx+1)*aR(nx+1)/(aR(nx+1)-aL(nx+1)))*(enR(nx+1)-enL(nx+1)))*dxi
      
     den = -den
   end subroutine CUMethod
   
   !> @brief Kurganov-Tadmor central scheme
   subroutine KTMethod(dro,dpx,den,roL,pxL,enL,roR,pxR,enR,x)
     real(wp),dimension(-2:nx+3) :: roL, pxL, enL, roR, pxR, enR, x
     real(wp), dimension(-2:nx+2) :: dro, dpx, den, A, pL, pR
     real(wp) :: p, cL, cR, dxi
     integer :: i
     
     dxi = 1.0_wp/dx

! Set spectral radius
     do i=1,nx+1
        pL(i)=gamma1*(enL(i) - half*(pxL(i)*pxL(i))/roL(i))
        pR(i)=gamma1*(enR(i) - half*(pxR(i)*pxR(i))/roR(i))
        if ((pL(i)<0).or.(pR(i)<0)) then 
           print*, "Error: negative pressure in KTMethod:",pL(i),pR(i)
           pL(i) = max(pL(i), 1e-10_wp)
           pR(i) = max(pR(i), 1e-10_wp)
        endif
     
        cL=SQRT(gamma*pL(i)/roL(i))
        cR=SQRT(gamma*pR(i)/roR(i))
        A(i)=MAX(ABS(pxL(1)/roL(1)-cL),ABS(pxR(1)/roR(1)-cR),&
               ABS(pxL(1)/roL(1)+cL),ABS(pxR(1)/roR(1)+cR))
     enddo
     
! flux of ro: ro*u
      dro(1)=-((pxL(1)) + (pxR(1)) - A(1)*(roR(1) - roL(1)))
      do i=2,nx
        p = ((pxL(i)) + (pxR(i)) - A(i)*(roR(i) - roL(i)))
        dro(i-1) = dro(i-1) + p
        dro(i  ) = dro(i  ) - p
      enddo
     dro(nx)=dro(nx)+((pxL(nx+1)) + (pxR(nx+1)) - A(nx+1)*(roR(nx+1) - roL(nx+1)))

     dro(:) = -dro(:)*dxi



! flux of ro*u: ro*u*u+p
    dpx(1)=-((pxL(1)**2/roL(1)+pL(1))+(pxR(1)**2/roR(1)+pR(1))&
         -A(1)*(pxR(1)-pxL(1)))
     do i=2,nx
        p = (pxL(i)**2/roL(i) + pL(i)) + (pxR(i)**2/roR(i) + pR(i)) - A(i)*(pxR(i) - pxL(i))
        dpx(i-1) = dpx(i-1) + p
        dpx(i  ) = dpx(i  ) - p
     enddo
     dpx(nx)=dpx(nx) + (pxL(nx+1)**2/roL(nx+1) + pL(nx+1)) +&
                           (pxR(nx+1)**2/roR(nx+1) + pR(nx+1)) - A(nx+1)*(pxR(nx+1) - pxL(nx+1))
    
	 dpx(:) = -dpx(:)*dxi

! flux ef energy: u*(E+p)
     den(1) = -((pxL(1)/roL(1)*(enL(1) + pL(1))) + (pxR(1)/roR(1)*(enR(1)+pR(1))) - A(1)*(enR(1) - enL(1)))
     do i=2,nx
        p = ((pxL(i)/roL(i)*(enL(i) + pL(i))) + (pxR(i)/roR(i)*(enR(i)+pR(i))) - A(i)*(enR(i) - enL(i)))
        den(i-1) = den(i-1) + p
        den(i  ) = den(i  ) - p
     enddo
     den(nx) = den(nx) + ((pxL(nx+1)/roL(nx+1)*(enL(nx+1) + pL(nx+1))) + (pxR(nx+1)/roR(nx+1)*(enR(nx+1)+pR(nx+1)))&
                 - A(nx+1)*(enR(nx+1) - enL(nx+1)))
	 
     den(:) = -den(:)*dxi
     
   end subroutine KTMethod

end program


