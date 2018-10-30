module stream_module

  use amrex_fort_module, only : amrex_real, dim=>amrex_spacedim
  
  implicit none

  public

contains

  subroutine pushvtog(lo, hi, dlo, dhi, U, U_lo, U_hi, nc) bind(c,name='pushvtog')
    implicit none
    integer, intent(in) :: nc, lo(3),  hi(3), dlo(3), dhi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    real(amrex_real), intent(inout) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3),nc)

    integer :: n
    real(amrex_real) :: xlo(3), dx(3)

    ! Make up something for these that gets what we want
    dx(1:3)  = 1._amrex_real
    xlo(1:3) = 0._amrex_real

    do n = 1,nc
       call hoextraptocc(U(:,:,:,n),U_lo(1),U_lo(2),U_lo(3),U_hi(1),U_hi(2),U_hi(3),lo,hi,dx,xlo)
    enddo
      
  end subroutine pushvtog

  subroutine vtrace(T, T_lo, T_hi, nT, loc, loc_lo, loc_hi, nl,&
     &     ids, n_ids, g, g_lo, g_hi, computeVec, strm, strm_lo, strm_hi,&
     &     ncs, dx, plo, hRK, errFlag) bind(C,name="vtrace")
    implicit none
    integer, intent(in) ::  nT, nl, computeVec, ncs, n_ids
    integer, intent(in) :: T_lo(3),T_hi(3),g_lo(3),g_hi(3),loc_lo(3),loc_hi(3),strm_lo(3),strm_hi(3)
    real(amrex_real), intent(in) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3),nT)
    real(amrex_real), intent(inout) :: g(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),dim)
    real(amrex_real), intent(in) :: loc(loc_lo(1):loc_hi(1),loc_lo(2):loc_hi(2),loc_lo(3):loc_hi(3),nl)
    real(amrex_real), intent(inout) :: strm(strm_lo(1):strm_hi(1),strm_lo(2):strm_hi(2),strm_lo(3):strm_hi(3),ncs)
    integer, intent(in) :: ids(0:n_ids-1)
    real(amrex_real), intent(in) :: dx(dim), plo(dim), hRK

    integer :: i,j,k,n,m
    real(amrex_real) :: x(3),xp(3),xm(3)

    logical :: ok
    integer :: myproc, errFlag

    call bl_pd_myproc(myproc)

    ! Compute gradient field assuming grow cells valid
    if (computeVec .eq. 1) then
       do k=g_lo(3),g_hi(3)
          do j=g_lo(2),g_hi(2)
             do i=g_lo(1),g_hi(1)
                  g(i,j,k,1) = T(i+1,j,k,1) - T(i-1,j,k,1)
                  g(i,j,k,2) = T(i,j+1,k,1) - T(i,j-1,k,1)
                  if (dim.eq.3) then
                     g(i,j,k,3) = T(i,j,k+1,1) - T(i,j,k-1,1)
                  endif
             enddo
          enddo
       enddo
    endif

    ! Loop over list of nodes in this box and do work
    !  loc contains all the nodes, ids is a list of those in this box
    !  strm holds path trace data: i=cnt of node in this box, j=index away from iso, k=unused, m=path comp
    !  Note that the ids come from the inside_nodes structure which is 1-based, but the
    !  indexing into strm is boxlib 0-based, subtract 1 when setting j here.
    errFlag = 0
    do i=0,n_ids-1
       j = ids(i) - 1
       x(1:dim) = loc(j,loc_lo(2),loc_lo(3),1:dim)
       strm(i,0,strm_lo(3),1:dim) = x(1:dim)
       do m=1,nT
          call ntrpv(x,T(T_lo(1),T_lo(2),T_lo(3),m),T_lo,T_hi,dx,plo,&
               strm(i,0,strm_lo(3),dim+m),1,ok)
          if (ok .eqv. .false.) then
             errFlag = 1
             return
          endif
       enddo

       ! Integrate with RK4, interpolate input data at new positions
       xp(1:dim) = x(1:dim)
       xm(1:dim) = x(1:dim)
       do n=-1,strm_lo(2),-1
          call RK4(xm,-hRK,g,g_lo,g_hi,dx,plo,ok)
          if (ok .eqv. .false.) then
             errFlag = 1
             return
          endif
          strm(i,n,strm_lo(3),1:dim) = xm(1:dim)
          do m=1,nT
             call ntrpv(xm,T(T_lo(1),T_lo(2),T_lo(3),m),T_lo,T_hi,dx,plo,&
                  strm(i,n,strm_lo(3),dim+m),1,ok)
             if (ok .eqv. .false.) then
                errFlag = 1
                return
             endif
          enddo
       enddo
         
       do n=1,strm_hi(2)
          call RK4(xp,+hRK,g,g_lo,g_hi,dx,plo,ok)
          if (ok .eqv. .false.) then
             errFlag = 1
             return
          endif
          strm(i,n,strm_lo(3),1:dim) = xp(1:dim)
          do m=1,nT
             call ntrpv(xp,T(T_lo(1),T_lo(2),T_lo(3),m),T_lo,T_hi,dx,plo,&
                  strm(i,n,strm_lo(3),dim+m),1,ok)
             if (ok .eqv. .false.) then
                errFlag = 1
                return
             endif
          enddo
       enddo
    enddo
  end subroutine vtrace

  subroutine RK4(x,h,g,g_lo,g_hi,dx,plo,ok)
    implicit none
    real(amrex_real), intent(inout) :: x(3)
    real(amrex_real), intent(in) :: h,dx(3),plo(3)
    integer, intent(in) :: g_lo(3),g_hi(3)
    real(amrex_real), intent(in) :: g(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),3)

    real(amrex_real) :: xx(3),vec(3),k1(3),k2(3),k3(3),k4(3)
    logical ok

    xx(1:dim) = x(1:dim)
    call ntrpv(xx,g,g_lo,g_hi,dx,plo,vec,dim,ok)
    if ( ok .eqv. .false. ) return
    call vnrml(vec)

    k1(1:dim) = vec(1:dim)*h
    xx(1:dim) = x(1:dim) + k1(1:dim)*.5d0
    call ntrpv(xx,g,g_lo,g_hi,dx,plo,vec,dim,ok)
    if ( ok .eqv. .false. ) return
    call vnrml(vec)
      
    k2(1:dim) = vec(1:dim)*h
    xx(1:dim) = x(1:dim) + k2(1:dim)*.5d0
    call ntrpv(xx,g,g_lo,g_hi,dx,plo,vec,dim,ok)
    if ( ok .eqv. .false. ) return
    call vnrml(vec)
    
    k3(1:dim) = vec(1:dim)*h
    xx(1:dim) = x(1:dim) + k3(1:dim)
    call ntrpv(xx,g,g_lo,g_hi,dx,plo,vec,dim,ok)
    if ( ok .eqv. .false. ) return
    call vnrml(vec)
     
    k4(1:dim) = vec(1:dim)*h
    x(1:dim) = x(1:dim) + (k1(1:dim)+k4(1:dim))/6.d0 + (k2(1:dim)+k3(1:dim))/3.d0
  end subroutine RK4

  subroutine ntrpv(x,g,g_lo,g_hi,dx,plo,u,nc,ok)
    implicit none
    integer, intent(in) :: nc
    real(amrex_real), intent(in) :: x(dim), dx(dim), plo(dim)
    real(amrex_real), intent(inout) :: u(nc)
    integer, intent(in) :: g_lo(dim), g_hi(dim)
    real(amrex_real), intent(in) :: g(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),nc)
    logical, intent(inout) :: ok

    integer :: b(dim), i
    real(amrex_real) :: n(dim), tmp
    do i=1,dim
       tmp = (x(i) - plo(i)) / dx(i) - 0.5d0
       b(i) = FLOOR( tmp )
       n(i) = ( x(i) - ( (b(i)+0.5d0)*dx(i) + plo(i) ) )/dx(i)
       n(i) = MAX(0.d0,MIN(1.d0,n(i)))
    enddo

    ok = .true.
    do i=1,dim
       if (b(i).lt.g_lo(i) .or. b(i).gt.g_hi(i)) ok = .false.
    enddo
    if (.not. ok) then
       print *,'b:',b
       print *,'DIMS:',g_lo,g_hi
       ok = .false.
       return
    endif

    do i=1,nc
       if (dim.eq.2) then
          u(i) = &
               +       n(1)  *    n(2)   * g(b(1)+1,b(2)+1,g_lo(3),i) &
               +       n(1)  *(1.d0-n(2))* g(b(1)+1,b(2)  ,g_lo(3),i) &
               +  (1.d0-n(1))*    n(2)   * g(b(1)  ,b(2)+1,g_lo(3),i) &
               +  (1.d0-n(1))*(1.d0-n(2))* g(b(1)  ,b(2)  ,g_lo(3),i)
       else
          u(i) = &
               +       n(1)  *    n(2)   *    n(3)    * g(b(1)+1,b(2)+1,b(3)+1,i) &
               +       n(1)  *(1.d0-n(2))*    n(3)    * g(b(1)+1,b(2)  ,b(3)+1,i) &
               +       n(1)  *    n(2)   *(1.d0-n(3)) * g(b(1)+1,b(2)+1,b(3)  ,i) &
               +       n(1)  *(1.d0-n(2))*(1.d0-n(3)) * g(b(1)+1,b(2)  ,b(3)  ,i) &
               +  (1.d0-n(1))*    n(2)   *    n(3)    * g(b(1)  ,b(2)+1,b(3)+1,i) &
               +  (1.d0-n(1))*(1.d0-n(2))*    n(3)    * g(b(1)  ,b(2)  ,b(3)+1,i) &
               +  (1.d0-n(1))*    n(2)   *(1.d0-n(3)) * g(b(1)  ,b(2)+1,b(3)  ,i) &
               +  (1.d0-n(1))*(1.d0-n(2))*(1.d0-n(3)) * g(b(1)  ,b(2)  ,b(3)  ,i)
       endif
    enddo
  end subroutine ntrpv

  subroutine vnrml(vec)
    implicit none
    real(amrex_real), intent(inout) :: vec(dim)
    real(amrex_real) :: eps, sum
    parameter (eps=1.e-12)
    integer :: i

    sum = 0._amrex_real
    do i=1,dim
       sum = sum + vec(i)*vec(i)
    enddo
    if (sum .gt. eps) then
       sum = sqrt(sum)
       do i=1,dim
          vec(i) = vec(i) / sum
       enddo
    endif
  end subroutine vnrml

end module stream_module
