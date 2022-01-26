module process_module

  use amrex_fort_module, only : amrex_real
  
  implicit none

  public

contains

  subroutine process(lo, hi, dlo, dhi, U, U_lo, U_hi, nc, plo, dx) bind(c,name='process')
    implicit none
    integer, intent(in) :: nc, lo(3),  hi(3), dlo(3), dhi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    real(amrex_real), intent(inout) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3),nc)
    real(amrex_real), intent(in) :: plo(3), dx(3)

  end subroutine process

end module process_module
