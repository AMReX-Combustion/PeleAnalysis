module util_module

  public :: get_spec_names, init_mech

contains

  subroutine init_mech() bind(C, name="init_mech")
    use network, only: network_init
    use transport_module, only: transport_init
    implicit none
    call network_init()
    call transport_init()
  end subroutine init_mech

  function num_species() bind(C, name="num_species") result(L)
    use network, only : nspecies
    implicit none
    integer :: L
    L = nspecies
  end function num_species
  
  function num_elements() bind(C, name="num_elements") result(L)
    use network, only : nelements
    implicit none
    integer :: L
    L = nelements
  end function num_elements
  
  function max_spec_namelen() bind(C, name="max_spec_namelen") result(L)
    use network, only : L_spec_name
    implicit none
    integer :: L
    L = L_spec_name
  end function max_spec_namelen
  
  function max_elem_namelen() bind(C, name="max_elem_namelen") result(L)
    use network, only : L_elem_name
    implicit none
    integer :: L
    L = L_elem_name
  end function max_elem_namelen
  
  subroutine get_spec_names(spec_names_out,ispec,len) bind(C, name="get_spec_names")
    use network, only : spec_names, L_spec_name
    implicit none
    integer, intent(in   ) :: ispec
    integer, intent(inout) :: len
    integer, intent(inout) :: spec_names_out(len)
    integer   :: i

    if (len .lt. L_spec_name) then
       print *,'Result container too short for spec name'
       call bl_pd_abort()
    endif

    len = len_trim(spec_names(ispec+1))

    do i = 1,len
       spec_names_out(i) = ichar(spec_names(ispec+1)(i:i))
    end do

  end subroutine get_spec_names

  subroutine get_elem_names(elem_names_out,ielem,len) bind(C, name="get_elem_names")
    use network, only : elem_names, L_elem_name
    implicit none
    integer, intent(in   ) :: ielem
    integer, intent(inout) :: len
    integer, intent(inout) :: elem_names_out(len)
    integer   :: i

    if (len .lt. L_elem_name) then
       print *,'Result container too short for elem name'
       call bl_pd_abort()
    endif

    len = len_trim(elem_names(ielem+1))

    do i = 1,len
       elem_names_out(i) = ichar(elem_names(ielem+1)(i:i))
    end do

  end subroutine get_elem_names

end module util_module
