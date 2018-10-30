module grad_module
  
  use amrex_fort_module, only : amrex_real, dim=>amrex_spacedim
  use amrex_bc_types_module
  use amrex_filcc_module
  
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

  subroutine gradient(lo, hi, U, U_lo, U_hi, G, G_lo, G_hi, dx) bind(c,name='gradient')
    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3), G_lo(3), G_hi(3)
    real(amrex_real), intent(in   ) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3))
    real(amrex_real), intent(inout) :: G(G_lo(1):G_hi(1),G_lo(2):G_hi(2),G_lo(3):G_hi(3),dim+1)
    real(amrex_real), intent(in)    :: dx(3)
    real(amrex_real) :: dxInv(3)
    integer :: i,j,k

    dxInv(1:3) = 0.5_amrex_real / dx(1:3)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             G(i,j,k,1) = dxInv(1)*(U(i+1,j,k) - U(i-1,j,k))
             G(i,j,k,2) = dxInv(2)*(U(i,j+1,k) - U(i,j-1,k))
             if (dim .gt. 2) then
                G(i,j,k,3) = dxInv(3)*(U(i,j,k+1) - U(i,j,k-1))
             else
                G(i,j,k,3) = 0._amrex_real
             endif
             G(i,j,k,dim+1) = SQRT(G(i,j,k,1)**2 + G(i,j,k,2)**2 + G(i,j,k,3)**2)
          enddo
       enddo
    enddo
  end subroutine gradient

end module grad_module
