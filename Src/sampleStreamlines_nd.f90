module stream_module

  use amrex_fort_module, only : amrex_real, dim=>amrex_spacedim
  
  implicit none

  public

contains

  subroutine interpstream(loc, loc_lo, loc_hi, nl,&
       fab, fab_lo, fab_hi, np, strm, strm_lo, strm_hi, dx, plo) bind(C,name="interpstream")
    
    implicit none
    integer, intent(in) ::  nl, np
    integer, intent(in) :: loc_lo(3),loc_hi(3),fab_lo(3),fab_hi(3),strm_lo(3),strm_hi(3)
    real(amrex_real), intent(in) :: loc(loc_lo(1):loc_hi(1),loc_lo(2):loc_hi(2),loc_lo(3):loc_hi(3),nl)
    real(amrex_real), intent(in) :: fab(fab_lo(1):fab_hi(1),fab_lo(2):fab_hi(2),fab_lo(3):fab_hi(3),np)
    real(amrex_real), intent(inout) :: strm(strm_lo(1):strm_hi(1),strm_lo(2):strm_hi(2),strm_lo(3):strm_hi(3),np)
    real(amrex_real), intent(in) :: dx(3), plo(3)

    integer :: i,j,k
    real(amrex_real) :: x(3),vtmp(3)
    logical :: ok

    do k=loc_lo(3),loc_hi(3)
       do j=0,loc_lo(2),-1
          do i=loc_lo(1),loc_hi(1)
             x(1:dim) = loc(i,j,k,1:dim)
             call ntrpv(x,fab,fab_lo,fab_hi,dx,plo,vtmp,np,ok)
             if (.not. ok) then
                if (j.eq.0) then
                   call bl_abort('Seed not in valid region for interp')
                else
                   call bl_abort('Interp bad, increase nGrow')
                endif
             else
                strm(i,j,k,1:np) = vtmp(1:np)
             endif
          enddo
       enddo
       do j=1,loc_hi(2)
          do i=loc_lo(1),loc_hi(1)
             x(1:dim) = loc(i,j,k,1:dim)
             call ntrpv(x,fab,fab_lo,fab_hi,dx,plo,vtmp,np,ok)
             if (.not. ok) then
                call bl_abort('Interp bad, increase nGrow')
             else
                strm(i,j,k,1:np) = vtmp(1:np)
             endif
          enddo
       enddo
    enddo
  end subroutine interpstream

  subroutine ntrpv(x,g,g_lo,g_hi,dx,plo,u,nc,ok)
    implicit none
    integer, intent(in) :: nc
    real(amrex_real), intent(in) :: x(3), dx(3), plo(3)
    real(amrex_real), intent(inout) :: u(nc)
    integer, intent(in) :: g_lo(3), g_hi(3)
    real(amrex_real), intent(in) :: g(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),nc)
    logical, intent(inout) :: ok

    integer :: b(3), i
    real(amrex_real) :: n(3), tmp
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

  subroutine set_distance(loc, loc_lo, loc_hi, res, res_lo, res_hi) bind(C, name='set_distance')
    implicit none
    integer, intent(in) :: loc_lo(3),loc_hi(3),res_lo(3),res_hi(3)
    real(amrex_real), intent(inout) :: loc(loc_lo(1):loc_hi(1),loc_lo(2):loc_hi(2),loc_lo(3):loc_hi(3),3)
    real(amrex_real), intent(inout) :: res(res_lo(1):res_hi(1),res_lo(2):res_hi(2),res_lo(3):res_hi(3))

    integer :: i,j,k
    real(amrex_real) :: d,dx,dy,dz

    do k=res_lo(3),res_hi(3)
       do i=res_lo(1),res_hi(1)

          d = 0.d0
          res(i,0,0) = d
          do j=-1,res_lo(2),-1
             dx = loc(i,j,k,1) - loc(i,j+1,k,1)
             dy = loc(i,j,k,2) - loc(i,j+1,k,2)
             if (dim.eq.2) then
                dz = 0._amrex_real
             else
                dz = loc(i,j,k,3) - loc(i,j+1,k,3)
             endif
             d = d + SQRT(dx*dx + dy*dy + dz*dz)
             res(i,j,k) = -d
          enddo

          d = 0.d0
          do j=1,res_hi(2)
             dx = loc(i,j,k,1) - loc(i,j-1,k,1)
             dy = loc(i,j,k,2) - loc(i,j-1,k,2)
             if (dim.eq.2) then
                dz = 0._amrex_real
             else
                dz = loc(i,j,k,3) - loc(i,j-1,k,3)
             endif
             d = d + SQRT(dx*dx + dy*dy + dz*dz)
             res(i,j,k) = d
          enddo
       enddo
    enddo
  end subroutine set_distance
end module stream_module
