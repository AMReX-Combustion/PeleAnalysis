module util_module

  public :: get_spec_names, init_mech

contains

  subroutine init_mech() bind(C, name="init_mech")

    use network, only: network_init
    use transport_module, only: transport_init

    call network_init()
    call transport_init()

  end subroutine init_mech

  subroutine get_spec_names(spec_names_out,ispec,len) &
       bind(C, name="get_spec_names")

    use network, only : spec_names

    implicit none

    integer, intent(in   ) :: ispec
    integer, intent(inout) :: len
    integer, intent(inout) :: spec_names_out(len)

    ! Local variables
    integer   :: i

    len = len_trim(spec_names(ispec+1))

    do i = 1,len
       spec_names_out(i) = ichar(spec_names(ispec+1)(i:i))
    end do

  end subroutine get_spec_names

end module util_module
