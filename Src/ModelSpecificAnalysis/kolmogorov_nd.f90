module grad_module
  
  use amrex_fort_module, only : amrex_real, dim=>amrex_spacedim
  use amrex_bc_types_module
  use amrex_filcc_module
  
  implicit none

  real(amrex_real), parameter :: viscosity=1.85e-1 ![cm2/s]
  integer :: id_lo(3) = (/0,0,0/)
  integer :: id_hi(3) = (/0,0,0/)

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

  subroutine compute_Kolmogorov(lo, hi, U, U_lo, U_hi, Um, Um_lo, Um_hi, eta, eta_lo, eta_hi, dx) bind(c,name='compute_Kolmogorov')
    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3), Um_lo(3), Um_hi(3), eta_lo(3), eta_hi(3) 
    real(amrex_real), intent(in   ) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3),3)
    real(amrex_real), intent(in   ) :: Um(Um_lo(1):Um_hi(1),Um_lo(2):Um_hi(2),Um_lo(3):Um_hi(3),3)
    real(amrex_real), intent(inout) :: eta(eta_lo(1):eta_hi(1),eta_lo(2):eta_hi(2),eta_lo(3):eta_hi(3),3)
    real(amrex_real), intent(in)    :: dx(3)
    real(amrex_real) :: dxInv(3)
    real(amrex_real), dimension(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) :: epsln
    real(amrex_real), dimension(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), 3, 3) :: grad_u, tau 
    real(amrex_real), dimension(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3),3) :: U_prime
    integer :: i,j,k
    integer :: nx, ny, nz
    real(amrex_real) :: mu_cubed = viscosity**3

    nx = hi(1)-lo(1)
    ny = hi(2)-lo(2)
    nz = hi(3)-lo(3)
    id_hi(:) = hi(:)
    id_lo(:) = lo(:)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             U_prime(i,j,k,1) = U(i,j,k,1) - Um(i,j,k,1)
             U_prime(i,j,k,2) = U(i,j,k,2) - Um(i,j,k,2)
             U_prime(i,j,k,3) = U(i,j,k,3) - Um(i,j,k,3)
          enddo
       enddo
    enddo 

    call calculate_dissipation(Um, Um_lo, Um_hi, U, U_lo, U_hi, dx, epsln, grad_u, tau, lo, hi)

    ! print*,"Computing grad for the first time"
    ! do i=1,3
    !    call gradient(Um(Um_lo(1):Um_hi(1),Um_lo(2):Um_hi(2),Um_lo(3):Um_hi(3),i), dx, grad_u(id_lo(1):id_hi(1),id_lo(2):id_hi(2),id_lo(3):id_hi(3),1:3, i))
    ! enddo
    ! print*,"Exited successfully!"
    ! tau = grad_u


    dxInv(1:3) = 0.5_amrex_real / dx(1:3)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ! epsln(i,j,k)       =   tau(i,j,k,1,1)*grad_u(i,j,k,1,1)    &
             !                      + tau(i,j,k,1,2)*grad_u(i,j,k,2,1)    &
             !                      + tau(i,j,k,1,3)*grad_u(i,j,k,3,1)    &
             !                      + tau(i,j,k,2,1)*grad_u(i,j,k,1,2)    &
             !                      + tau(i,j,k,2,2)*grad_u(i,j,k,2,2)    &
             !                      + tau(i,j,k,2,3)*grad_u(i,j,k,3,2)    &
             !                      + tau(i,j,k,3,1)*grad_u(i,j,k,1,3)    &
             !                      + tau(i,j,k,3,2)*grad_u(i,j,k,2,3)    &
             !                      + tau(i,j,k,3,3)*grad_u(i,j,k,3,3)    
             eta(i,j,k,1) = (mu_cubed/epsln(i,j,k))**0.25
             eta(i,j,k,2) = epsln(i,j,k)
             eta(i,j,k,3) = grad_u(i,j,k,1,1)
          enddo
       enddo
    enddo    

  end subroutine compute_Kolmogorov

  subroutine gradient(U, U_lo, U_hi, dx, grad, grad_lo, grad_hi)
    implicit none
    ! integer, intent(in) :: nx, ny, nz
    integer, intent(in) :: U_lo(3),  U_hi(3), grad_lo(3),  grad_hi(3)
    real(amrex_real), intent(in   ) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3))
    real(amrex_real), intent(inout) :: grad(grad_lo(1):grad_hi(1),grad_lo(2):grad_hi(2),grad_lo(3):grad_hi(3), 3)
    real(amrex_real), intent(in)    :: dx(3)
    real(amrex_real) :: dxInv(3)
    integer :: i,j,k
    ! integer :: nx, ny, nz

    ! nx = hi(1)-lo(1)
    ! ny = hi(2)-lo(2)
    ! nz = hi(3)-lo(3)

    dxInv(1:3) = 0.5_amrex_real / dx(1:3)
    do k = id_lo(3),id_hi(3)
       do j = id_lo(2),id_hi(2)
          do i = id_lo(1),id_hi(1)
             grad(i,j,k,1) = dxInv(1)*(U(i+1,j,k) - U(i-1,j,k))
             grad(i,j,k,2) = dxInv(2)*(U(i,j+1,k) - U(i,j-1,k))
             grad(i,j,k,3) = dxInv(3)*(U(i,j,k+1) - U(i,j,k-1))
             ! G(i,j,k,dim+1) = SQRT(G(i,j,k,1)**2 + G(i,j,k,2)**2 + G(i,j,k,3)**2)
          enddo
       enddo
    enddo
  end subroutine gradient

subroutine computeStressTensor(grad_u,lo,hi)

!----------------------------------------------------------------------
!Rewritten by Ramanan Sankaran - 01/05/05
!  To save memory this function was rewritten.
!  Now Tau is returned overwriting grad_u
!  Just so the code is readable, i have a pointer tau=>grad_u
!  But remember that , now, both are the same memory locations
!  That is why we need to momentarily store as grad_u_tmp
!----------------------------------------------------------------------
!
 !--------------------------------------------------------------------------------
 ! given the velocity gradient tensor, this routine computes the stress tensor
 ! for a newtonian fluid
 !
 ! In index notation, the stress tensor is written as:
 !
 !  tau_ij = mu [ d(u_i)/dx_j + d(u_j)/dx_i ]
 !         - del_ij [ 2/3(mu-K) * d(u_k)/dx_k ]
 !
 ! where:
 !    mu  =  viscosity
 !    K   =  bulk viscosity (ASSUMED TO BE ZERO HERE - true for monatomic gases)
 !    del_ij = Kronicker delta (1 if i=j, else 0)
 !--------------------------------------------------------------------------------

 ! use param_m,  only: nx,ny,nz
integer, intent(in) :: lo(3),  hi(3)
 real(amrex_real), intent(inout),  target, dimension(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3,3) :: grad_u
 real(amrex_real), pointer, dimension(:,:,:,:,:) :: tau

 real(amrex_real) :: grad_u_tmp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3,3)
 real(amrex_real) :: sumterm(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

 integer :: i,j,k,m,n

 tau => grad_u


 grad_u_tmp = grad_u

 sumterm = 0.0
 do n = 1, 3
   sumterm(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = sumterm(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) + grad_u_tmp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),n,n)
 end do

 ROW: do m=1,3
    tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),m,m) = 2.0*viscosity * (grad_u_tmp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),m,m) - sumterm(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))/3.0)

    ! copmute the rest of the columns in this row, using
    ! symmetry for elements to the left of the diagonal.
    ! really we don't need to store the whole thing, but 
    ! this makes life much easier...
    COLUMN: do n=m+1,3
       tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),m,n) = viscosity * ( grad_u_tmp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),m,n) + grad_u_tmp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),n,m) )
       tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),n,m) = tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),m,n)
    enddo COLUMN
 enddo ROW


 return
end subroutine computeStressTensor

subroutine calculate_dissipation(big_u, big_u_lo, big_u_hi, lil_u, lil_u_lo, lil_u_hi, dx, epsln, grad_u, tau, lo, hi)

!----------------------------------------
! - routine calculates the dissipation of kinetic energy
! - Needs the big_u which is mean+prime velocities for calculating
!   the stress tensor.
! - Needs little u to calculate the gradient tensor, which is then 
!   contracted against the stress tensor.
!   Little u is big_u - favre_mean_u
!----------------------------------------
! use transport_m, only: computestresstensor
implicit none

integer, intent(in) :: lo(3),  hi(3)
integer, intent(in) :: big_u_lo(3),  big_u_hi(3), lil_u_lo(3),  lil_u_hi(3)
real(amrex_real), dimension(big_u_lo(1):big_u_hi(1),big_u_lo(2):big_u_hi(2),big_u_lo(3):big_u_hi(3), 3), intent(in) :: big_u
real(amrex_real), dimension(lil_u_lo(1):lil_u_hi(1),lil_u_lo(2):lil_u_hi(2),lil_u_lo(3):lil_u_hi(3), 3), intent(in) :: lil_u
real(amrex_real), intent(in)    :: dx(3)
! real, dimension(nx, ny, nz), intent(in) :: volum
real(amrex_real), dimension(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), intent(inout) :: epsln
real(amrex_real), dimension(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), 3, 3), intent(inout) :: grad_u, tau
integer :: i

! Compute a grad_u based on the big_u
! call computevectorgradient(big_u, grad_u)
 do i=1,3
    call gradient(big_u(big_u_lo(1):big_u_hi(1),big_u_lo(2):big_u_hi(2),big_u_lo(3):big_u_hi(3),i), &
     big_u_lo(1:3), big_u_hi(1:3), dx, grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:3,i),lo(1:3),hi(1:3))
 enddo

tau = grad_u
call computeStressTensor(tau,lo,hi)

! Recompute grad_u based on the lil_u
! call computevectorgradient(lil_u, grad_u)
 do i=1,3
    call gradient(lil_u(lil_u_lo(1):lil_u_hi(1),lil_u_lo(2):lil_u_hi(2),lil_u_lo(3):lil_u_hi(3),i), &
     lil_u_lo(1:3), lil_u_hi(1:3), dx, grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:3,i),lo(1:3),hi(1:3))
 enddo

! contract the stress with the velocity gradients to form the dissipation
! epsilon=tau_{ij} du_j/dx_i

  epsln(:,:,:)       =   tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1,1)*grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1,1)    &
                       + tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1,2)*grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2,1)    &
                       + tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1,3)*grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3,1)    &
                       + tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2,1)*grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1,2)    &
                       + tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2,2)*grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2,2)    &
                       + tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2,3)*grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3,2)    &
                       + tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3,1)*grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1,3)    &
                       + tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3,2)*grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2,3)    &
                       + tau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3,3)*grad_u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3,3)    

! divide by density
  ! epsln(:,:,:)=epsln(:,:,:)

return
end subroutine calculate_dissipation
end module grad_module
