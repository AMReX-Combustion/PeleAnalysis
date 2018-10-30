module isoval_module

  use amrex_fort_module, only : amrex_real, dim=>amrex_spacedim
  
  implicit none

  public

contains

  subroutine setloc(lo, hi, U, U_lo, U_hi, dx, xlo) bind(c,name="setloc")
    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    real(amrex_real), intent(in) :: xlo(3), dx(3)
    real(amrex_real), intent(inout) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3),1:dim)

    integer :: i,j,k
    real(amrex_real) :: loc(3)

    do k=lo(3),hi(3)
       loc(3) = (k + 0.5d0)*dx(3) + xlo(3)
       do j=lo(2),hi(2)
          loc(2) = (j + 0.5d0)*dx(2) + xlo(2)
          do i=lo(1),hi(1)
             loc(1) = (i + 0.5d0)*dx(1) + xlo(1)
             U(i,j,k,1:dim) = loc(1:dim)
          enddo
       enddo
    enddo

  end subroutine setloc

  subroutine setcloc(lo, hi, U, U_lo, U_hi, fdx, xlo, ratio) bind(c,name="setcloc")
    implicit none
    integer, intent(in) :: lo(3),  hi(3), ratio
    integer, intent(in) :: U_lo(3), U_hi(3)
    real(amrex_real), intent(in) :: xlo(3), fdx(3)
    real(amrex_real), intent(inout) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3),1:dim)

    integer :: i,j,k,ci,cj,ck
    real(amrex_real) :: loc(3)

    do k=lo(3),hi(3)
       ck = k/ratio
       if (k.lt.0) ck = ck-1
       loc(3) = (ck + 0.5d0)*fdx(3)*ratio + xlo(3)
       do j=lo(2),hi(2)
          cj = j/ratio
          if (j.lt.0) cj = cj-1
          loc(2) = (cj + 0.5d0)*fdx(2)*ratio + xlo(2)
          do i=lo(1),hi(1)
             ci = i/ratio
             if (i.lt.0) ci = ci-1
             loc(1) = (ci + 0.5d0)*fdx(1)*ratio + xlo(1)
             U(i,j,k,1:dim) = loc(1:dim)
          enddo
       enddo
    enddo

  end subroutine setcloc

end module isoval_module
