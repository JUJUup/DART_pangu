
! DART note: this file started life as:
! /glade/p/cesm/cseg/collections/cesm1_1_1/models/atm/cam/src/dynamics/homme/dyn_comp.F90

Module dyn_comp

  use shr_kind_mod, only: r8 => shr_kind_r8
  use domain_mod, only : domain1d_t
  use element_mod, only : element_t, elem_state_t
  use time_mod, only : TimeLevel_t
  use hybvcoord_mod, only : hvcoord_t
  use hybrid_mod, only : hybrid_t
  use perf_mod, only: t_startf, t_stopf
  use cam_logfile, only : iulog
  use time_manager,     only: is_first_step
  use spmd_utils,  only : iam, npes_cam => npes
  use pio,         only: file_desc_t
  use cslam_control_volume_mod, only : cslam_struct

  implicit none
  private


  ! PUBLIC MEMBER FUNCTIONS:
  public dyn_init1, dyn_init2, dyn_run, dyn_final

  ! PUBLIC DATA MEMBERS:
  public dyn_import_t, dyn_export_t


  type (TimeLevel_t)   , public :: TimeLevel     ! main time level struct (used by tracers)

!  type (elem_state_t), save, target :: dyn_state

  type dyn_import_t
     type (element_t), pointer :: elem(:)
     type (cslam_struct), pointer :: cslam(:)
  end type dyn_import_t

  type dyn_export_t
     type (element_t), pointer :: elem(:)
     type (cslam_struct), pointer :: cslam(:)
  end type dyn_export_t
  type (hvcoord_t) :: hvcoord
  integer, parameter  ::  DYN_RUN_SUCCESS           = 0
  integer, parameter  ::  DYN_RUN_FAILURE           = -1

  ! !DESCRIPTION: This module implements the HOMME Dynamical Core as
  !               an ESMF gridded component.  It is specific to HOMME
  !               and does not use ESMF.
  !
  ! \paragraph{Overview}
  !
  !   This module contains an ESMF wrapper for the Homme
  !   Dynamical Core used in the Community Atmospheric Model. 
  !
  ! !REVISION HISTORY:
  !
  !  JPE  06.05.31:  created
  !
  !----------------------------------------------------------------------

  ! Enumeration of DYNAMICS_IN_COUPLINGS


  logical, parameter         :: DEBUG = .true.

  real(r8), parameter        :: ONE    = 1.0D0

  character(*), parameter, public :: MODULE_NAME = "dyn_comp"
  character(*), parameter, public :: VERSION     = "$Id$" 
  type (domain1d_t), pointer, public :: dom_mt(:)      


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dyn_init1(fh, NLFileName, dyn_in, dyn_out)

  ! Initialize the dynamical core

    use pio,              only: file_desc_t
    use hycoef,           only: hycoef_init
    use ref_pres,         only: ref_pres_init

    use pmgrid,           only: dyndecomp_set
    use dyn_grid,         only: dyn_grid_init, cslam, elem, set_horiz_grid_cnt_d, &
                                get_dyn_grid_parm, set_horiz_grid_cnt_d
    use rgrid,            only: fullgrid
    use spmd_utils,       only: mpi_integer, mpicom, mpi_logical
    use spmd_dyn,         only: spmd_readnl
    use interpolate_mod,  only: interpolate_analysis
    use native_mapping,   only: create_native_mapping_files, native_mapping_readnl

    use dimensions_mod,   only: globaluniquecols, nelem, nelemd, nelemdmax
    use prim_driver_mod,  only: prim_init1
    use thread_mod,       only: nthreads
    use parallel_mod,     only: par, initmp
    use namelist_mod,     only: readnl
    use control_mod,      only: runtype

    ! PARAMETERS:
    type(file_desc_t),   intent(in)  :: fh       ! PIO file handle for initial or restart file
    character(len=*),    intent(in)  :: NLFileName
    type (dyn_import_t), intent(OUT) :: dyn_in
    type (dyn_export_t), intent(OUT) :: dyn_out

#ifdef _OPENMP    
    integer omp_get_num_threads
#endif
    integer :: neltmp(3)
    logical :: nellogtmp(7)
    integer :: npes_homme

    ! Initialize dynamics grid
    call dyn_grid_init()

    ! Read in the number of tasks to be assigned to Homme (needed by initmp)
    call spmd_readnl(NLFileName, npes_homme)
    ! Initialize the homme structure that holds the MPI decomposition information
    par=initmp(npes_homme)

    ! Read the homme specific part of the namelist
    call readnl(par, NLFileName)

    ! override the setting in the homme namelist, it's redundent anyway
    if (.not. is_first_step()) runtype = 1

    ! Initialize hybrid coordinate arrays.
    call hycoef_init(fh)

    ! Initialize physics grid reference pressures (needed by initialize_radbuffer)
    call ref_pres_init()

    ! legacy reduced grid code -- should be removed
    fullgrid=.true.

#ifdef _OPENMP    
!   Set by driver
!$omp parallel
    nthreads = omp_get_num_threads()
!$omp end parallel
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: number of OpenMP threads = ", nthreads
       write(iulog,*) " "
    endif
#if defined (ELEMENT_OPENMP)
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: using OpenMP within element instead of across elements"
       write(iulog,*) " "
    endif
#endif
#else
    nthreads = 1
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: openmp not activated"
       write(iulog,*) " "
    endif
#endif
    if(iam < par%nprocs) then
       call prim_init1(elem,cslam,par,dom_mt,TimeLevel)

       dyn_in%elem => elem
       dyn_out%elem => elem
       dyn_in%cslam => cslam
       dyn_out%cslam => cslam
    
       call set_horiz_grid_cnt_d(GlobalUniqueCols)


       neltmp(1) = nelemdmax
       neltmp(2) = nelem
       neltmp(3) = get_dyn_grid_parm('plon')
       nellogtmp(1:7) = interpolate_analysis(1:7)
    else
       nelemd = 0
       neltmp(1) = 0
       neltmp(2) = 0
       neltmp(3) = 0
       nellogtmp(1:7) = .true.
    endif

    dyndecomp_set = .true.



    if (par%nprocs .lt. npes_cam) then
! Broadcast quantities to auxiliary processes
       call mpibcast(neltmp, 3, mpi_integer, 0, mpicom)
       call mpibcast(nellogtmp, 7, mpi_logical, 0, mpicom)
       if (iam .ge. par%nprocs) then
          nelemdmax = neltmp(1)
          nelem     = neltmp(2)
          call set_horiz_grid_cnt_d(neltmp(3))
          interpolate_analysis(1:7) = nellogtmp(1:7)
       endif
    endif


    !
    ! This subroutine creates mapping files using homme basis functions if requested
    !
    call native_mapping_readnl(NLFileName)
    call create_native_mapping_files( par, elem,'native')
    call create_native_mapping_files( par, elem,'bilin')

  end subroutine dyn_init1


  subroutine dyn_init2(dyn_in)
    use dimensions_mod,   only: nlev, nelemd
    use prim_driver_mod,  only: prim_init2, prim_run
    use prim_si_ref_mod,  only: prim_set_mass
    use hybrid_mod,       only: hybrid_create
    use hycoef,           only: hyam, hybm, hyai, hybi, ps0
    use parallel_mod,     only: par
    use time_manager,     only: dtime,get_nstep   ! physics timestep
    use time_mod,         only: se_nsplit=>nsplit, tstep, time_at
    use control_mod,      only: moisture, runtype, qsplit
    use thread_mod,       only: nthreads
    use thread_mod,       only: omp_get_thread_num
    use cam_control_mod,  only: aqua_planet, ideal_phys, adiabatic
    use comsrf,           only: landm, sgh, sgh30
    use nctopo_util_mod,  only: nctopo_util_driver
! KDR to prevent multiple HommeMapping writes
    use cam_instance,     only: inst_index


    type (dyn_import_t), intent(inout) :: dyn_in

    type(element_t),    pointer :: elem(:)
    type(cslam_struct), pointer :: cslam(:)

    integer :: ithr, nets, nete, ie, k
    real(r8), parameter :: Tinit=300.0_r8
    real(r8) :: dyn_ps0
    type(hybrid_t) :: hybrid

    !
    !  Note: dtime = progress made in one timestep.  value in namelist
    !        dtime = the frequency at which physics is called
    !        tstep = the dynamics timestep:  
    !
    !    Leapfrog looks like:   u(3) = u(1) + 2*tstep*u(2)   
    !    u(1) = time-tstep 
    !    u(2) = time
    !    u(3) = time+tstep 
    !  
    !    Physics looks like:    u(3) = u(1) + dt_phys*PHYSICS(U(1))
    !  
    ! so with se_nsplit=1:    dtime=tstep   dt_phys=2*tstep
    !
    ! In general:  dtime=se_nsplit*tstep,  dt_phys=se_nsplit*tstep + tstep
    ! 

    elem  => dyn_in%elem
    cslam => dyn_in%cslam

    dyn_ps0=ps0/100.D0
    hvcoord%hyam=hyam
    hvcoord%hyai=hyai
    hvcoord%hybm=hybm
    hvcoord%hybi=hybi
    hvcoord%ps0=dyn_ps0  
    do k=1,nlev
       hvcoord%hybd(k) = hvcoord%hybi(k+1) - hvcoord%hybi(k)
    end do
    if(iam < par%nprocs) then

#if (! defined ELEMENT_OPENMP)
       !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ie,ithr,nets,nete,hybrid)
#endif
       ithr=omp_get_thread_num()
       nets=dom_mt(ithr)%start
       nete=dom_mt(ithr)%end
       hybrid = hybrid_create(par,ithr,NThreads)

       tstep = dtime/float(se_nsplit*qsplit)

       moisture='moist'

       if(adiabatic) then
          moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%q(:,:,:,:,:)=0.0D0
                elem(ie)%derived%fq(:,:,:,:,:)=0.0D0
             end do
          end if
       else if(ideal_phys) then
          moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%lnps(:,:,:) =LOG(dyn_ps0)

                elem(ie)%state%ps_v(:,:,:) =dyn_ps0

                elem(ie)%state%phis(:,:)=0.0D0

                elem(ie)%state%T(:,:,:,:) =Tinit

                elem(ie)%state%v(:,:,:,:,:) =0.0D0

                elem(ie)%state%q(:,:,:,:,:)=0.0D0

             end do
          end if
       else if(aqua_planet .and. runtype==0)  then
          do ie=nets,nete
             !          elem(ie)%state%lnps(:,:,:) =LOG(dyn_ps0)
             !          elem(ie)%state%ps_v(:,:,:) =dyn_ps0
             elem(ie)%state%phis(:,:)=0.0D0
          end do
          if(allocated(landm)) landm=0.0_r8
          if(allocated(sgh)) sgh=0.0_r8
          if(allocated(sgh30)) sgh30=0.0_r8
       end if

       do ie=nets,nete
          elem(ie)%derived%FM=0.0D0
          elem(ie)%derived%FT=0.0D0
          elem(ie)%derived%FQ=0.0D0
       end do

       ! initial homme (subcycled) nstep
       TimeLevel%nstep = get_nstep()*se_nsplit*qsplit

       ! scale PS to achieve prescribed dry mass
       if (runtype == 1) then
          ! exact restart
          TimeLevel%nstep0=TimeLevel%nstep+1
       else
          ! new run, scale mass to value given in namelist, if needed
          call prim_set_mass(elem, TimeLevel,hybrid,hvcoord,nets,nete)  
          TimeLevel%nstep0=2   ! This will be the first full leapfrog step
       endif
       call prim_init2(elem,cslam,hybrid,nets,nete, TimeLevel, hvcoord)
       !
       ! This subroutine is used to create nc_topo files, if requested
       ! 
       call nctopo_util_driver(elem,hybrid,nets,nete)
#if (! defined ELEMENT_OPENMP)
       !$OMP END PARALLEL 
#endif
    end if

! KDR; there should be a restriction that only instance X write out this file.
    if(inst_index==1) then
       call write_grid_mapping(par, elem)
    end if
    
  end subroutine dyn_init2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  RUN --- Driver for the 
  !
  ! !INTERFACE:
  subroutine dyn_run( dyn_state, rc )

    ! !USES:
    use parallel_mod,     only : par
    use prim_driver_mod,  only: prim_run, prim_run_subcycle
    use dimensions_mod,   only : nlev
    use thread_mod,       only: omp_get_thread_num, nthreads
    use time_mod,         only: se_nsplit=>nsplit, tstep
    use control_mod,      only: tstep_type
    use hybrid_mod,       only: hybrid_create
!    use perf_mod, only : t_startf, t_stopf
    implicit none


    type (dyn_export_t), intent(inout)       :: dyn_state   !  container
    type(hybrid_t) :: hybrid

    integer, intent(out)               :: rc      ! Return code
    integer ::  n
    integer :: nets, nete, ithr
    integer :: ie
    real(r8) :: tstep_tmp

    ! !DESCRIPTION:
    !
    if(iam < par%nprocs) then
#if (! defined ELEMENT_OPENMP)
       !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid,n)
#endif
       ithr=omp_get_thread_num()
       nets=dom_mt(ithr)%start
       nete=dom_mt(ithr)%end
       hybrid = hybrid_create(par,ithr,NThreads)

       do n=1,se_nsplit
          if (tstep_type==1) then
             ! forward-in-time RK, with subcycling
             call prim_run_subcycle(dyn_state%elem,dyn_state%cslam,hybrid,nets,nete,&
                  tstep, TimeLevel, hvcoord)
          else
             ! leapfrog
             call prim_run(dyn_state%elem, hybrid,nets,nete, tstep, TimeLevel, hvcoord, "leapfrog")
          endif
       end do


#if (! defined ELEMENT_OPENMP)
       !$OMP END PARALLEL
#endif
    end if
    rc = DYN_RUN_SUCCESS

    !EOC
  end subroutine dyn_run
  !-----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dyn_final(DYN_STATE, RESTART_FILE)

    type (elem_state_t), target     :: DYN_STATE
    character(LEN=*)   , intent(IN) :: RESTART_FILE



  end subroutine dyn_final



  subroutine write_grid_mapping(par, elem)
    use parallel_mod,     only: parallel_t
    use element_mod, only : element_t
    use cam_pio_utils, only : cam_pio_createfile, pio_subsystem
    use pio, only : file_desc_t, pio_def_dim, var_desc_t, pio_int, pio_def_var, &
         pio_enddef, pio_closefile, pio_initdecomp, io_desc_t, pio_write_darray, &
         pio_freedecomp, pio_setdebuglevel
    use dimensions_mod, only : np, nelem, nelemd
    use dof_mod, only : createmetadata

    type(parallel_t) :: par
    type(element_t) :: elem(:)
    type(file_desc_t) :: nc
    type(var_desc_t) :: vid
    type(io_desc_t) :: iodesc
    integer :: dim1, dim2, ierr, i, j, ie, cc, base, ii, jj
    integer, parameter :: npm12 = (np-1)*(np-1)
    integer :: subelement_corners(npm12*nelemd,4)
    integer :: dof(npm12*nelemd*4)


    ! Create a CS grid mapping file for postprocessing tools

       ! write meta data for physics on GLL nodes
       call cam_pio_createfile(nc, 'HommeMapping.nc', 0)
   
       ierr = pio_def_dim(nc, 'ncenters', npm12*nelem, dim1)
       ierr = pio_def_dim(nc, 'ncorners', 4, dim2)
       ierr = pio_def_var(nc, 'element_corners', PIO_INT, (/dim1,dim2/),vid)
    
       ierr = pio_enddef(nc)
       call createmetadata(par, elem, subelement_corners)

       jj=0
       do cc=0,3
          do ie=1,nelemd
             base = ((elem(ie)%globalid-1)+cc*nelem)*npm12
             ii=0
             do j=1,np-1
                do i=1,np-1
                   ii=ii+1
                   jj=jj+1
                   dof(jj) = base+ii
                end do
             end do
          end do
       end do

       call pio_initdecomp(pio_subsystem, pio_int, (/nelem*npm12,4/), dof, iodesc)

       call pio_write_darray(nc, vid, iodesc, reshape(subelement_corners,(/nelemd*npm12*4/)), ierr)
       
       call pio_freedecomp(nc, iodesc)
       
       call pio_closefile(nc)

  end subroutine write_grid_mapping

end module dyn_comp



