module util_module

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
  
  function num_reactions() bind(C, name="num_reactions") result(L)
    use network, only : nreactions
    implicit none
    integer :: L
    L = nreactions
  end function num_reactions
  
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

  subroutine get_rns_with_X_on_L(rns,len,Xid) bind(C, name="get_rns_with_X_on_L")
    use network, only : nreactions, nspecies
    use fuego_chemistry
    implicit none
    integer, intent(in   ) :: Xid
    integer, intent(inout) :: len
    integer, intent(inout) :: rns(len)
    integer :: max_num_spec, num_spec, ireac, lenout, ispec
    integer, allocatable :: specIDs(:), specNUs(:)

    call CKINU(0,max_num_spec,specIDs,specNUs)
    allocate(specIDs(1:max_num_spec))
    allocate(specNUs(1:max_num_spec))

    lenout = 0
    do ireac = 1,nreactions
       call CKINU(ireac,num_spec,specIDs,specNUs)
       if (num_spec .lt. 0  .or.  num_spec .gt. nspecies) call bl_pd_abort()
       do ispec = 1,num_spec
          if (specIDs(ispec).eq.Xid .and. specNUs(ispec).lt.0) then
             if (lenout.eq.len) call bl_pd_abort() ! Not enough space to write output
             lenout = lenout + 1
             rns(lenout) = ireac
           endif
       enddo
    enddo
    len = lenout
     deallocate(specIDs,specNUs)
    
  end subroutine get_rns_with_X_on_L

  subroutine get_rns_with_X_on_R(rns,len,Xid) bind(C, name="get_rns_with_X_on_R")
    use network, only : nreactions, nspecies
    use fuego_chemistry
    implicit none
    integer, intent(in   ) :: Xid
    integer, intent(inout) :: len
    integer, intent(inout) :: rns(len)
    integer :: max_num_spec, num_spec, ireac, lenout, ispec
    integer, allocatable :: specIDs(:), specNUs(:)

    call CKINU(0,max_num_spec,specIDs,specNUs)
    allocate(specIDs(1:max_num_spec))
    allocate(specNUs(1:max_num_spec))

    lenout = 0
    do ireac = 1,nreactions
       call CKINU(ireac,num_spec,specIDs,specNUs)
       if (num_spec .lt. 0  .or.  num_spec .gt. nspecies) call bl_pd_abort()
       do ispec = 1,num_spec
          if (specIDs(ispec).eq.Xid .and. specNUs(ispec).gt.0) then
             if (lenout.eq.len) call bl_pd_abort() ! Not enough space to write output
             lenout = lenout + 1
             rns(lenout) = ireac
           endif
       enddo
    enddo
    len = lenout
     deallocate(specIDs,specNUs)
    
  end subroutine get_rns_with_X_on_R

  subroutine get_spec_nu_for_rn(ids,nu,len,ireac) bind(C, name="get_spec_nu_for_rn")
    use network, only : nreactions, nspecies
    use fuego_chemistry
    implicit none
    integer, intent(in   ) :: ireac
    integer, intent(inout) :: len
    integer, intent(inout) :: ids(len), nu(len)
    integer :: max_num_spec, num_spec, ispec
    integer, allocatable :: specIDs(:), specNUs(:)

    call CKINU(0,max_num_spec,specIDs,specNUs)
    allocate(specIDs(1:max_num_spec))
    allocate(specNUs(1:max_num_spec))

    call CKINU(ireac,num_spec,specIDs,specNUs)
    do ispec=1,num_spec
       ids(ispec) = specIDs(ispec)
       nu(ispec) = specNUs(ispec)
    enddo
    len = num_spec
    deallocate(specIDs,specNUs)    
  end subroutine get_spec_nu_for_rn

  subroutine get_rns_map(rmap) bind(C, name="get_rns_map")
    use network, only : nreactions
    use fuego_chemistry
    implicit none
    integer, intent(inout) :: rmap(nreactions)
    integer :: ireac

    call GET_REACTION_MAP(rmap)
    
  end subroutine get_rns_map

  subroutine get_elt_comp(ispec,comp) bind(C, name="get_elt_comp")
    use network, only : nspecies, nelements
    use fuego_chemistry
    implicit none
    integer, intent(in) :: ispec
    integer, intent(inout) :: comp(*)
    integer :: ncf(nspecies*nelements)
    integer :: i
    call CKNCF(nelements,ncf)
    do i=1,nelements
       comp(i) = ncf((ispec-1)*nelements + i)
    enddo
  end subroutine get_elt_comp

end module util_module
