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

  subroutine normalize(lo, hi, U, U_lo, U_hi, S, S_lo, S_hi, nmin, nmax, dx) bind(C,name='normalize')
    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3), S_lo(3), S_hi(3)
    real(amrex_real), intent(in   ) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3))
    real(amrex_real), intent(inout) :: S(S_lo(1):S_hi(1),S_lo(2):S_hi(2),S_lo(3):S_hi(3))
    real(amrex_real) :: nmin, nmax, dx(3), fac
    real(amrex_real) :: x, y, z
    integer :: i,j, k
    fac = 1._amrex_real / (nmax - nmin)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             S(i,j,k) = ( U(i,j,k) - nmin ) * fac
             ! uncomment to replace T with a field of r for testing
             !x = (i+0.5)*dx(1)
             !y = (j+0.5)*dx(2)
             !z = (k+0.5)*dx(3)
             !S(i,j,k) = SQRT((x-.015)**2 + (y-.015)**2 + (z-.015)**2 )
          end do
       end do
    end do
  end subroutine normalize

  subroutine mcurv(lo, hi, T, T_lo, T_hi, curv, curv_lo, curv_hi,&
       wrk, wrk_lo, wrk_hi, delta, sym) bind(C,name='mcurv')
    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: T_lo(3),T_hi(3),curv_lo(3),curv_hi(3),wrk_lo(3),wrk_hi(3)
    real(amrex_real), intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
    real(amrex_real), intent(inout) :: curv(curv_lo(1):curv_hi(1),curv_lo(2):curv_hi(2),curv_lo(3):curv_hi(3))
    real(amrex_real), intent(inout) :: wrk(wrk_lo(1):wrk_hi(1),wrk_lo(2):wrk_hi(2),wrk_lo(3):wrk_hi(3),dim)
    real(amrex_real), intent(inout) :: delta(dim), sym(dim)

    integer :: i,j,k
    real(amrex_real) :: mag,gT(dim),dxInv(dim),halfDxInv(dim),dxInvSq(dim)
    real(amrex_real) :: Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Txz,eps,zero
    parameter (eps = 1.d-8)
    parameter (zero = 0._amrex_real)

    dxInv(1:dim) = 1._amrex_real / delta(1:dim)
    dxInvSq(1:dim) = dxInv(1:dim) * dxInv(1:dim)
    halfDxInv(1:dim) = 0.5_amrex_real * dxInv(1:dim)

    if (dim.eq.2) then
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)+1
                gT(1) = (T(i,j,k) + T(i,j-1,k) - T(i-1,j,k) - T(i-1,j-1,k)) * halfDxInv(1)
                gT(2) = (T(i,j,k) - T(i,j-1,k) + T(i-1,j,k) - T(i-1,j-1,k)) * halfDxInv(2)
                mag = SQRT(gT(1)*gT(1) + gT(2)*gT(2))

                if (mag.eq.0._amrex_real) then
                   wrk(i,j,k,1:dim) = 0._amrex_real
                else
                   wrk(i,j,k,1:dim) = -gT(1:dim)/mag ! Normal points into fuel
                endif
             enddo
          enddo
       enddo
       ! Get curvature at centers from nodel normal, get normal at centers
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                curv(i,j,k) = .5_amrex_real*((wrk(i+1,j+1,k,1)-wrk(i,j+1,k,1)) * dxInv(1) &
                     +                       (wrk(i+1,j  ,k,1)-wrk(i,j  ,k,1)) * dxInv(1) &
                     +                       (wrk(i+1,j+1,k,2)-wrk(i+1,j,k,2)) * dxInv(2) &
                     +                       (wrk(i  ,j+1,k,2)-wrk(i  ,j,k,2)) * dxInv(2))
             enddo
          enddo
       enddo
    else      
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                Tx = halfDxInv(1)*(T(i+1,j,k) - T(i-1,j,k))
                Ty = halfDxInv(2)*(T(i,j+1,k) - T(i,j-1,k))
                Tz = halfDxInv(3)*(T(i,j,k+1) - T(i,j,k-1))
                Txx = dxInvSq(1)*(T(i+1,j,k) - 2*T(i,j,k) + T(i-1,j,k))
                Tyy = dxInvSq(2)*(T(i,j+1,k) - 2*T(i,j,k) + T(i,j-1,k))
                Tzz = dxInvSq(3)*(T(i,j,k+1) - 2*T(i,j,k) + T(i,j,k-1))
                Txy = halfDxInv(1)*halfDxInv(2)*(T(i+1,j+1,k) - T(i-1,j+1,k) - T(i+1,j-1,k) + T(i-1,j-1,k))
                Txz = halfDxInv(1)*halfDxInv(3)*(T(i+1,j,k+1) - T(i-1,j,k+1) - T(i+1,j,k-1) + T(i-1,j,k-1))
                Tyz = halfDxInv(2)*halfDxInv(3)*(T(i,j+1,k+1) - T(i,j-1,k+1) - T(i,j+1,k-1) + T(i,j-1,k-1))
                if (sym(1).eq.1) then
                   Tx  = zero
                   Txx = zero
                   Txy = zero
                   Txz = zero
                endif
                if (sym(2).eq.1) then
                   Ty  = zero
                   Tyy = zero
                   Txy = zero
                   Tyz = zero
                endif
                if (sym(3).eq.1) then
                   Tz  = zero
                   Tzz = zero
                   Txz = zero
                   Tyz = zero
                endif

                mag = max(1.0d-12, SQRT(Tx*Tx + Ty*Ty + Tz*Tz))

                curv(i,j,k) = -0.5_amrex_real*(Txx + Tyy + Tzz &
                     - ( Tx*(Tx*Txx + Ty*Txy + Tz*Txz) &
                     +   Ty*(Tx*Txy + Ty*Tyy + Tz*Tyz) &
                     +   Tz*(Tx*Txz + Ty*Tyz + Tz*Tzz) )/mag**2 )/mag

                ! uncomment to get back average radius of curvature instead of mean curvature
                !curv(i,j,k) = sign( max(abs(curv(i,j,k)),1.e-12), curv(i,j,k))
                !curv(i,j,k) = -1.d0 / curv(i,j,k)
                
             enddo
          enddo
       enddo
    endif
  end subroutine mcurv

  subroutine strainrate(lo, hi, u, u_lo, u_hi, T, T_lo, T_hi, sr, sr_lo, sr_hi,&
     &     wrk, wrk_lo, wrk_hi, delta) bind(C, name='strainrate')

    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: u_lo(3),u_hi(3),T_lo(3),T_hi(3),sr_lo(3),sr_hi(3),wrk_lo(3),wrk_hi(3)
    real(amrex_real), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),dim)
    real(amrex_real), intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
    real(amrex_real), intent(inout) :: sr(sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3))
    real(amrex_real), intent(inout) :: wrk(wrk_lo(1):wrk_hi(1),wrk_lo(2):wrk_hi(2),wrk_lo(3):wrk_hi(3),dim)
    real(amrex_real), intent(in   ) :: delta(dim)

    integer :: i,j,k,ii,jj,kk,c
    real(amrex_real) :: dxInv(dim), halfDxInv(dim)
    real(amrex_real) :: mag,gT(dim),eps,ddi2,ddj2,ddi3,ddj3,ddk3,str2,str3

    ddi2(ii,jj,kk,c)= halfDxInv(1)*(u(ii,jj,kk,c)-u(ii-1,jj,kk,c)+u(ii,jj-1,kk,c)-u(ii-1,jj-1,kk,c))
    ddj2(ii,jj,kk,c)= halfDxInv(2)*(u(ii,jj,kk,c)-u(ii,jj-1,kk,c)+u(ii-1,jj,kk,c)-u(ii-1,jj-1,kk,c))
    str2(ii,jj,kk) &
         = ddi2(ii,jj,kk,1)*wrk(ii,jj,kk,1)*wrk(ii,jj,kk,1) &
         + ddi2(ii,jj,kk,2)*wrk(ii,jj,kk,1)*wrk(ii,jj,kk,2) &
         + ddj2(ii,jj,kk,1)*wrk(ii,jj,kk,2)*wrk(ii,jj,kk,1) &
         + ddj2(ii,jj,kk,2)*wrk(ii,jj,kk,2)*wrk(ii,jj,kk,2)
    ddi3(ii,jj,kk,c)&
         = 0.25*(u(ii,jj,kk  ,c)-u(ii-1,jj,kk  ,c)+u(ii,jj-1,kk  ,c)-u(ii-1,jj-1,kk  ,c) &
         +       u(ii,jj,kk-1,c)-u(ii-1,jj,kk-1,c)+u(ii,jj-1,kk-1,c)-u(ii-1,jj-1,kk-1,c)) * dxInv(1)
    ddj3(ii,jj,kk,c)&
         = 0.25*(u(ii,jj,kk  ,c)-u(ii,jj-1,kk  ,c)+u(ii-1,jj,kk  ,c)-u(ii-1,jj-1,kk  ,c) &
         +       u(ii,jj,kk-1,c)-u(ii,jj-1,kk-1,c)+u(ii-1,jj,kk-1,c)-u(ii-1,jj-1,kk-1,c)) * dxInv(2)
    ddk3(ii,jj,kk,c)&
         = 0.25*(u(ii,jj  ,kk,c)-u(ii,jj  ,kk-1,c)+u(ii-1,jj  ,kk,c)-u(ii-1,jj  ,kk-1,c) &
         +       u(ii,jj-1,kk,c)-u(ii,jj-1,kk-1,c)+u(ii-1,jj-1,kk,c)-u(ii-1,jj-1,kk-1,c)) * dxInv(3)

    str3(ii,jj,kk) &
         = -ddi3(ii,jj,kk,1)*wrk(ii,jj,kk,1)*wrk(ii,jj,kk,1) &
         -  ddi3(ii,jj,kk,2)*wrk(ii,jj,kk,1)*wrk(ii,jj,kk,2) &
         -  ddi3(ii,jj,kk,3)*wrk(ii,jj,kk,1)*wrk(ii,jj,kk,3) &
         -  ddj3(ii,jj,kk,1)*wrk(ii,jj,kk,2)*wrk(ii,jj,kk,1) &
         -  ddj3(ii,jj,kk,2)*wrk(ii,jj,kk,2)*wrk(ii,jj,kk,2) &
         -  ddj3(ii,jj,kk,3)*wrk(ii,jj,kk,2)*wrk(ii,jj,kk,3) &
         -  ddk3(ii,jj,kk,1)*wrk(ii,jj,kk,3)*wrk(ii,jj,kk,1) &
         -  ddk3(ii,jj,kk,2)*wrk(ii,jj,kk,3)*wrk(ii,jj,kk,2) &
         -  ddk3(ii,jj,kk,3)*wrk(ii,jj,kk,3)*wrk(ii,jj,kk,3) &
         +  ddi3(ii,jj,kk,1) &
         +  ddj3(ii,jj,kk,2) &
         +  ddk3(ii,jj,kk,3)

    eps = 1.d-20
    dxInv(1:dim) = 1._amrex_real / delta(1:dim)
    halfDxInv(1:dim) = 0.5_amrex_real * dxInv(1:dim)

    if (dim.eq.2) then
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)+1
             gT(1) = 0.5*(T(i,j,k)-T(i-1,j,k)+T(i,j-1,k)-T(i-1,j-1,k))/delta(1)
             gT(2) = 0.5*(T(i,j,k)-T(i,j-1,k)+T(i-1,j,k)-T(i-1,j-1,k))/delta(2)
             mag = MAX(eps,SQRT(gT(1)*gT(1) + gT(2)*gT(2)))
             wrk(i,j,k,1) = -gT(2)/mag
             wrk(i,j,k,2) =  gT(1)/mag
          enddo
       enddo

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sr(i,j,k) = 0.25d0*(str2(i,j,k)+str2(i+1,j,k)+str2(i,j+1,k)+str2(i+1,j+1,k))
          enddo
       enddo
    else
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)+1
                gT(1) = 0.25 * &
                     ( T(i,j,k-1)-T(i-1,j,k-1)+T(i,j-1,k-1)-T(i-1,j-1,k-1) &
                     + T(i,j,  k)-T(i-1,j,  k)+T(i,j-1,  k)-T(i-1,j-1,  k) ) * dxInv(1)
                gT(2) = 0.25 * &
                     ( T(i,j,k-1)-T(i,j-1,k-1)+T(i-1,j,k-1)-T(i-1,j-1,k-1) &
                     + T(i,j,  k)-T(i,j-1,  k)+T(i-1,j,  k)-T(i-1,j-1,  k) ) * dxInv(2)
                gT(3) = 0.25 * &
                     ( T(i,j,  k)+T(i-1,j,  k)+T(i,j-1,  k)+T(i-1,j-1,  k) &
                     - T(i,j,k-1)-T(i-1,j,k-1)-T(i,j-1,k-1)-T(i-1,j-1,k-1) ) * dxInv(3)
                mag = MAX(eps,SQRT(gT(1) *gT(1) + gT(2)*gT(2) + gT(3)*gT(3)))                
                if (mag.eq.eps) then
                   wrk(i,j,k,1:dim) = 0._amrex_real
                else
                   wrk(i,j,k,1:dim) = -gT(1:dim)/mag
                endif
             enddo
          enddo
       enddo
       ! Get strain rate, -nn:Grad(U)+Div(U)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                ! Get average value at each of 8 surrounding nodes
                sr(i,j,k) = 0.125_amrex_real &
                     * ( str3(i,j,k  )+str3(i+1,j,k  )+str3(i,j+1,k  )+str3(i+1,j+1,k  ) &
                     +   str3(i,j,k+1)+str3(i+1,j,k+1)+str3(i,j+1,k+1)+str3(i+1,j+1,k+1) )
             enddo
          enddo
       enddo
    endif
  end subroutine strainrate

  subroutine straintensor(lo, hi, u, u_lo, u_hi, sr, sr_lo, sr_hi, delta) bind(C, name='straintensor')
    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: u_lo(3),u_hi(3),sr_lo(3),sr_hi(3)
    real(amrex_real), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),dim)
    real(amrex_real), intent(inout) :: sr(sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3),9)
    real(amrex_real), intent(in   ) :: delta(dim)

    real(amrex_real) :: ddi2, ddj2, ddi3, ddj3, ddk3, halfDxInv(dim), fourthDxInv(dim)
    integer :: i,j,k,ii,jj,kk,c

    ddi2(ii,jj,kk,c)= halfDxInv(1)*(u(ii,jj,kk,c)-u(ii-1,jj,kk,c)+u(ii,jj-1,kk,c)-u(ii-1,jj-1,kk,c))
    ddj2(ii,jj,kk,c)= halfDxInv(2)*(u(ii,jj,kk,c)-u(ii,jj-1,kk,c)+u(ii-1,jj,kk,c)-u(ii-1,jj-1,kk,c))

    ddi3(ii,jj,kk,c)= fourthDxInv(1)*(u(ii,jj,kk  ,c)-u(ii-1,jj,kk  ,c)+u(ii,jj-1,kk  ,c)-u(ii-1,jj-1,kk  ,c) &
         +                            u(ii,jj,kk-1,c)-u(ii-1,jj,kk-1,c)+u(ii,jj-1,kk-1,c)-u(ii-1,jj-1,kk-1,c))
    ddj3(ii,jj,kk,c)= fourthDxInv(2)*(u(ii,jj,kk  ,c)-u(ii,jj-1,kk  ,c)+u(ii-1,jj,kk  ,c)-u(ii-1,jj-1,kk  ,c) &
         +                            u(ii,jj,kk-1,c)-u(ii,jj-1,kk-1,c)+u(ii-1,jj,kk-1,c)-u(ii-1,jj-1,kk-1,c))
    ddk3(ii,jj,kk,c)= fourthDxInv(3)*(u(ii,jj  ,kk,c)-u(ii,jj  ,kk-1,c)+u(ii-1,jj  ,kk,c)-u(ii-1,jj  ,kk-1,c) &
         +                            u(ii,jj-1,kk,c)-u(ii,jj-1,kk-1,c)+u(ii-1,jj-1,kk,c)-u(ii-1,jj-1,kk-1,c))

    halfDxInv(1:dim) = .5_amrex_real / delta(1:dim)
    fourthDxInv(1:dim) = .5_amrex_real * halfDxInv(1:dim)

    if (dim.eq.2) then
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                sr(i,j,k,1:9) = 0.d0
                sr(i,j,k,1) = 0.25d0*(ddi2(i,j,k,1)+ddi2(i+1,j,k,1)+ddi2(i,j+1,k,1)+ddi2(i+1,j+1,k,1))
                sr(i,j,k,2) = 0.25d0*(ddi2(i,j,k,2)+ddi2(i+1,j,k,2)+ddi2(i,j+1,k,2)+ddi2(i+1,j+1,k,2))
                sr(i,j,k,4) = 0.25d0*(ddj2(i,j,k,1)+ddj2(i+1,j,k,1)+ddj2(i,j+1,k,1)+ddj2(i+1,j+1,k,1))
                sr(i,j,k,5) = 0.25d0*(ddj2(i,j,k,2)+ddj2(i+1,j,k,2)+ddj2(i,j+1,k,2)+ddj2(i+1,j+1,k,2))
             enddo
          enddo
       enddo
    else
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
               sr(i,j,k,1) = 0.125d0* &
                    ( ddi3(i,j,k  ,1)+ddi3(i+1,j,k  ,1)+ddi3(i,j+1,k  ,1)+ddi3(i+1,j+1,k  ,1) &
                    + ddi3(i,j,k+1,1)+ddi3(i+1,j,k+1,1)+ddi3(i,j+1,k+1,1)+ddi3(i+1,j+1,k+1,1) )
               sr(i,j,k,2) = 0.125d0* &
                    ( ddi3(i,j,k  ,2)+ddi3(i+1,j,k  ,2)+ddi3(i,j+1,k  ,2)+ddi3(i+1,j+1,k  ,2) &
                    + ddi3(i,j,k+1,2)+ddi3(i+1,j,k+1,2)+ddi3(i,j+1,k+1,2)+ddi3(i+1,j+1,k+1,2) )
               sr(i,j,k,3) = 0.125d0* &
                    ( ddi3(i,j,k  ,3)+ddi3(i+1,j,k  ,3)+ddi3(i,j+1,k  ,3)+ddi3(i+1,j+1,k  ,3) &
                    + ddi3(i,j,k+1,3)+ddi3(i+1,j,k+1,3)+ddi3(i,j+1,k+1,3)+ddi3(i+1,j+1,k+1,3) )
               sr(i,j,k,4) = 0.125d0* &
                    ( ddj3(i,j,k  ,1)+ddj3(i+1,j,k  ,1)+ddj3(i,j+1,k  ,1)+ddj3(i+1,j+1,k  ,1) &
                    + ddj3(i,j,k+1,1)+ddj3(i+1,j,k+1,1)+ddj3(i,j+1,k+1,1)+ddj3(i+1,j+1,k+1,1) )
               sr(i,j,k,5) = 0.125d0* &
                    ( ddj3(i,j,k  ,2)+ddj3(i+1,j,k  ,2)+ddj3(i,j+1,k  ,2)+ddj3(i+1,j+1,k  ,2) &
                    + ddj3(i,j,k+1,2)+ddj3(i+1,j,k+1,2)+ddj3(i,j+1,k+1,2)+ddj3(i+1,j+1,k+1,2) )
               sr(i,j,k,6) = 0.125d0* &
                    ( ddj3(i,j,k  ,3)+ddj3(i+1,j,k  ,3)+ddj3(i,j+1,k  ,3)+ddj3(i+1,j+1,k  ,3) &
                    + ddj3(i,j,k+1,3)+ddj3(i+1,j,k+1,3)+ddj3(i,j+1,k+1,3)+ddj3(i+1,j+1,k+1,3) )
               sr(i,j,k,7) = 0.125d0* &
                    ( ddj3(i,j,k  ,1)+ddj3(i+1,j,k  ,1)+ddk3(i,j+1,k  ,1)+ddk3(i+1,j+1,k  ,1) &
                    + ddj3(i,j,k+1,1)+ddj3(i+1,j,k+1,1)+ddk3(i,j+1,k+1,1)+ddk3(i+1,j+1,k+1,1) )
               sr(i,j,k,8) = 0.125d0* &
                    ( ddj3(i,j,k  ,2)+ddj3(i+1,j,k  ,2)+ddk3(i,j+1,k  ,2)+ddk3(i+1,j+1,k  ,2) &
                    + ddj3(i,j,k+1,2)+ddj3(i+1,j,k+1,2)+ddk3(i,j+1,k+1,2)+ddk3(i+1,j+1,k+1,2) )
               sr(i,j,k,9) = 0.125d0* &
                    ( ddj3(i,j,k  ,3)+ddj3(i+1,j,k  ,3)+ddk3(i,j+1,k  ,3)+ddk3(i+1,j+1,k  ,3) &
                    + ddj3(i,j,k+1,3)+ddj3(i+1,j,k+1,3)+ddk3(i,j+1,k+1,3)+ddk3(i+1,j+1,k+1,3) )
            enddo
         enddo
      enddo
   endif
 end subroutine straintensor

 subroutine progressgrad(lo, hi, c, c_lo, c_hi, g, g_lo, g_hi, delta) bind(C, name='progressgrad')
    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: c_lo(3),c_hi(3),g_lo(3),g_hi(3)
    real(amrex_real), intent(in   ) :: c(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(amrex_real), intent(inout) :: g(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),dim)
    real(amrex_real), intent(in   ) :: delta(dim)

    real(amrex_real) :: ddi2,ddj2,ddi3,ddj3,ddk3
    integer :: i,j,k,ii,jj,kk,halfDxInv(dim),fourthDxInv(dim)

    ddi2(ii,jj,kk)= halfDxInv(1)*(c(ii,jj,kk)-c(ii-1,jj,kk)+c(ii,jj-1,kk)-c(ii-1,jj-1,kk))
    ddj2(ii,jj,kk)= halfDxInv(2)*(c(ii,jj,kk)-c(ii,jj-1,kk)+c(ii-1,jj,kk)-c(ii-1,jj-1,kk))
    ddi3(ii,jj,kk)= fourthDxInv(1)*(c(ii,jj,kk  )-c(ii-1,jj,kk  )+c(ii,jj-1,kk  )-c(ii-1,jj-1,kk  ) &
         +                          c(ii,jj,kk-1)-c(ii-1,jj,kk-1)+c(ii,jj-1,kk-1)-c(ii-1,jj-1,kk-1))
    ddj3(ii,jj,kk)= fourthDxInv(2)*(c(ii,jj,kk  )-c(ii,jj-1,kk  )+c(ii-1,jj,kk  )-c(ii-1,jj-1,kk  ) &
         +                          c(ii,jj,kk-1)-c(ii,jj-1,kk-1)+c(ii-1,jj,kk-1)-c(ii-1,jj-1,kk-1))
    ddk3(ii,jj,kk)= fourthDxInv(3)*(c(ii,jj  ,kk)-c(ii,jj  ,kk-1)+c(ii-1,jj  ,kk)-c(ii-1,jj  ,kk-1) &
         +                          c(ii,jj-1,kk)-c(ii,jj-1,kk-1)+c(ii-1,jj-1,kk)-c(ii-1,jj-1,kk-1))


    halfDxInv(1:dim) = .5_amrex_real / delta(1:dim)
    fourthDxInv(1:dim) = 0.5_amrex_real * halfDxInv(1:dim)

    if (dim.eq.2) then
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                g(i,j,k,1) = 0.25d0*(ddi2(i,j,k)+ddi2(i+1,j,k)+ddi2(i,j+1,k)+ddi2(i+1,j+1,k))
                g(i,j,k,2) = 0.25d0*(ddj2(i,j,k)+ddj2(i+1,j,k)+ddj2(i,j+1,k)+ddj2(i+1,j+1,k))
             enddo
          enddo
       enddo
    else
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                g(i,j,k,1) = 0.125d0* &
                     ( ddi3(i,j,k  )+ddi3(i+1,j,k  )+ddi3(i,j+1,k  )+ddi3(i+1,j+1,k  ) &
                     + ddi3(i,j,k+1)+ddi3(i+1,j,k+1)+ddi3(i,j+1,k+1)+ddi3(i+1,j+1,k+1) )
                g(i,j,k,2) = 0.125d0* &
                     ( ddj3(i,j,k  )+ddj3(i+1,j,k  )+ddj3(i,j+1,k  )+ddj3(i+1,j+1,k  ) &
                     + ddj3(i,j,k+1)+ddj3(i+1,j,k+1)+ddj3(i,j+1,k+1)+ddj3(i+1,j+1,k+1) )
                g(i,j,k,3) = 0.125d0* &
                     ( ddk3(i,j,k  )+ddk3(i+1,j,k  )+ddk3(i,j+1,k  )+ddk3(i+1,j+1,k  ) &
                     + ddk3(i,j,k+1)+ddk3(i+1,j,k+1)+ddk3(i,j+1,k+1)+ddk3(i+1,j+1,k+1) )
             enddo
          enddo
       enddo
    endif
  end subroutine progressgrad

  subroutine smooth(lo, hi, Tin, Tin_lo, Tin_hi, Tout, Tout_lo, Tout_hi) bind(C, name='smooth')
    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: Tin_lo(3),Tin_hi(3),Tout_lo(3),Tout_hi(3)
    real(amrex_real), intent(in   ) :: Tin(Tin_lo(1):Tin_hi(1),Tin_lo(2):Tin_hi(2),Tin_lo(3):Tin_hi(3))
    real(amrex_real), intent(inout) :: Tout(Tout_lo(1):Tout_hi(1),Tout_lo(2):Tout_hi(2),Tout_lo(3):Tout_hi(3))

    integer :: i,j,k,ii,jj,kk
    real(amrex_real) :: sixteenth, sixtyfourth
    parameter (sixteenth = 1._amrex_real / 16._amrex_real)
    parameter (sixtyfourth = 1._amrex_real / 64._amrex_real)

    if (dim.eq.2) then 
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                Tout(i,j,k) = 0._amrex_real
                do jj=0,1
                   do ii=0,1
                      Tout(i,j,k) = Tout(i,j,k) &
                           + Tin(i+ii,j+jj  ,k) + Tin(i+ii-1,j+jj  ,k) &
                           + Tin(i+ii,j+jj-1,k) + Tin(i+ii-1,j+jj-1,k)
                   enddo
                enddo
                Tout(i,j,k) = Tout(i,j,k) * sixteenth
             enddo
          enddo
       enddo
    else
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                Tout(i,j,k) = 0._amrex_real
                do kk=0,1
                   do jj=0,1
                      do ii=0,1
                         Tout(i,j,k) = Tout(i,j,k) &
                              + Tin(i+ii,j+jj,k+kk-1)   + Tin(i+ii-1,j+jj,k+kk-1) &
                              + Tin(i+ii,j+jj-1,k+kk-1) + Tin(i+ii-1,j+jj-1,k+kk-1) &
                              + Tin(i+ii,j+jj,k+kk)     + Tin(i+ii-1,j+jj,k+kk  ) &
                              + Tin(i+ii,j+jj-1,k+kk  ) + Tin(i+ii-1,j+jj-1,k+kk  )
                      enddo
                   enddo
                enddo
                Tout(i,j,k) = Tout(i,j,k) * sixtyfourth
             enddo
          enddo
       enddo
    endif
  end subroutine smooth

  subroutine gcurv(lo, hi, T, T_lo, T_hi, curv, curv_lo, curv_hi, delta) &
       bind(C, name='gcurv')
    implicit none
    integer, intent(in) :: lo(3),  hi(3)
    integer, intent(in) :: T_lo(3),T_hi(3),curv_lo(3),curv_hi(3)
    real(amrex_real), intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
    real(amrex_real), intent(inout) :: curv(curv_lo(1):curv_hi(1),curv_lo(2):curv_hi(2),curv_lo(3):curv_hi(3))
    real(amrex_real), intent(in) :: delta(dim)

    integer :: i,j,k
    real(amrex_real) :: Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Txz,eps,tdxI(3),dxSqI(3), A,B
    parameter (eps = 1.d-8)

    tdxI(1:dim) =  1._amrex_real / (2._amrex_real * delta(1:dim))
    dxSqI(1:dim) = 1._amrex_real / (delta(1:dim)*delta(1:dim))

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             Tx = tdxI(1)*(T(i+1,j,k) - T(i-1,j,k))
             Ty = tdxI(2)*(T(i,j+1,k) - T(i,j-1,k))
             Tz = tdxI(3)*(T(i,j,k+1) - T(i,j,k-1))

             Txx = dxSqI(1)*(T(i+1,j,k) - 2*T(i,j,k) + T(i-1,j,k))
             Tyy = dxSqI(2)*(T(i,j+1,k) - 2*T(i,j,k) + T(i,j-1,k))
             Tzz = dxSqI(3)*(T(i,j,k+1) - 2*T(i,j,k) + T(i,j,k-1))

             Txy = tdxI(1)*tdxI(2)*(T(i+1,j+1,k) - T(i-1,j+1,k) - T(i+1,j-1,k) + T(i-1,j-1,k))
             Txz = tdxI(1)*tdxI(3)*(T(i+1,j,k+1) - T(i-1,j,k+1) - T(i+1,j,k-1) + T(i-1,j,k-1))
             Tyz = tdxI(2)*tdxI(3)*(T(i,j+1,k+1) - T(i,j-1,k+1) - T(i,j+1,k-1) + T(i,j-1,k-1))

             A = 2._amrex_real*(Tx*Ty*(Txy*Tzz-Txz*Tyz) + Tx*Tz*(Txz*Tyy-Txy*Tyz) + Ty*Tz*(Tyz*Txx-Txy*Txz)) &
                  -         Tx*Tx*(Tyy*Tzz - Tyz*Tyz) - Ty*Ty*(Txx*Tzz-Txz*Txz) - Tz*Tz*(Txx*Tyy-Txy*Txy)
             B = Tx*Tx + Ty*Ty + Tz*Tz

             if (abs(B) .lt. eps) then
                curv(i,j,k) = 0._amrex_real
             else
                curv(i,j,k) = -A / B**2
             endif

          enddo
       enddo
    enddo
  end subroutine gcurv

end module stream_module
