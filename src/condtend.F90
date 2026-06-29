module condtend

   use phys_control, only: phys_getopts
   use chem_mods,    only: gas_pcnst
   use mo_tracname,  only: solsym
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_history,  only: outfld
   use physconst,    only: rair, gravit, pi
   use chem_mods,    only: adv_mass !molecular weights from mozart
   use ppgrid,       only: pcols, pver, pverp

   use aerosol_properties_mod,           only: aerosol_properties
   use sectional_aerosol_properties_mod, only: sectional_aerosol_properties
   use aerosol_state_mod,                only: aerosol_state, ptr2d_t
   use sectional_aerosol_state_mod,      only: sectional_aerosol_state

   implicit none
   private

!soa
! TODO: get rid of aero_sectional: from aero_props
! aero_props as input
! TODO: make everything, EVERYTHING!! allocatable!!
!   use aero_sectional,     only: secConstIndex

! Stuff to get rid of!
! secConstIndex(nspecies, nbins) is index for aerosols in qarray
! chemistryindex
! everything with "modal"
! ind_sec
! l_h2so4_chem = chemistryindex of h2so4

    real(r8), allocatable :: cond_sink_norm(:)       ![m3/#/s] condensation sink per particle in bin i
    real(r8), allocatable :: bin_centers(:)          ![m] bin centers
    ! [-] array of transformation of life cycle tracers
    integer :: l_h2so4_chem

contains

   subroutine registerCondensation()
      ! lifeCycleReceiver index setting for organics
      return
   end subroutine registerCondensation

!===============================================================================

   subroutine condensation_init(aero_props)

      !condensation coefficients:
      !Theory: Poling et al, "The properties of gases and liquids"
      !5th edition, eqn 11-4-4

      use cam_history,     only: addfld, add_default, fieldname_len, horiz_only

      ! dummy arguments
      type(sectional_aerosol_properties), intent(in) :: aero_props

      ! local
      real(r8), parameter :: aunit   = 1.6606e-27_r8  ! [kg] Atomic mass unit
      real(r8), parameter :: boltz   = 1.3806e-23_r8  ! [J/K/molec]
      real(r8), parameter :: t0      = 273.15_r8      ! [K] standard temperature
      real(r8), parameter :: p0      = 101325.0_r8    ! [Pa] Standard pressure
      real(r8), parameter :: radair  = 1.73e-10_r8    ! [m] Typical air molecule collision radius
      real(r8), parameter :: Mair    = 28.97_r8       ! [amu/molec] Molecular weight for dry air

      !Diffusion volumes for simple molecules [Poling et al], table 11-1
    ! TODO: vad seems like some kind of property for h2so4?
      real(r8), parameter :: vad     = 51.96_r8 ![cm3/mol]
      real(r8), parameter :: vadAir  = 19.7_r8                                          ![cm3/mol]
      real(r8), parameter :: aThird  = 1.0_r8/3.0_r8
      real(r8), parameter :: cm2Tom2 = 1.e-4_r8       ! convert from cm2 ==> m2

      !smb++ sectional
      real(r8), allocatable :: diff_coeff(:)   ! [m2/s] Diffusion coefficient sectional
      !smb-- sectional

      character(len=fieldname_len+3) :: fieldname_donor
      character(len=fieldname_len+3) :: fieldname_receiver
      character(128)                 :: long_name
      character(8)                   :: unit

      integer  :: iChem             !counter for chemical species
      integer  :: tracerIndex       !counter for chem. spec

      logical  :: history_aerosol
      logical  :: isAlreadyOnList(gas_pcnst)

      real(r8) :: mfv    ![m] mean free path
      real(r8) :: diff   ![m2/s] diffusion coefficient for cond. vap
      real(r8) :: molecularWeight !amu/molec molecular weight
      real(r8) :: Mdual  ![molec/amu] 1/M_1 + 1/M_2
      real(r8) :: rho    ![kg/m3] density of component in question
      real(r8) :: radmol ![m] radius molecule
      real(r8) :: th     !thermal velocity
      !smb++ sectional
      integer  :: ibin, ispecprop ! indexes
      integer  :: imozart
      integer  :: l_h2so4, l_h2so4_chem
      integer  :: sulfate_specprop_ndx ! index for sulfate species in aero_props
      character(len=16) :: type ! type of species in aero_props
      !smb-- sectional

      !-----------------------------------------------------------------------------------
      allocate(cond_sink_norm(aero_props%nbins()))
      allocate(diff_coeff(aero_props%nbins()))
      allocate(bin_centers(aero_props%nbins()))
      ! Couple the condenseable vapours to chemical species for properties and indexes
      ! add dimension for several species

      ! gas phase h2so4
      call cnst_get_ind('H2SO4',l_h2so4, abort=.true.)

      ! find first Mozart tracer imozart
      call cnst_get_ind(trim(solsym(1)), imozart, abort=.true.)

      ! Find the chemistry index for gas H2SO4
      l_h2so4_chem = l_h2so4 - imozart + 1

      do ispecprop = 1, aero_props%nspecies_tot() !TODO: make nicer to have several species
        type = aero_props%spectype(ispecprop)
        if (trim(type) == 'sulfate') then
            sulfate_specprop_ndx = ispecprop
            exit
        end if
      end do

      bin_centers = aero_props%bin_centers(aero_props%nbins())

      ! pick up densities and weights from aerosol properties
      rho    = aero_props%density(sulfate_specprop_ndx)                 ! aerosol type density
      molecularWeight = aero_props%molecular_weight(sulfate_specprop_ndx)    ! molecular weight of aerosol
      radmol = (3.0_r8*molecularWeight*aunit/(4.0_r8*pi*rho))**aThird        ! Radius of molecul (straight forward assuming spherical)
      Mdual  = 2.0_r8/(1.0_r8/Mair+1.0_r8/molecularWeight)                   ! factor of [1/m_1 + 1_m2]

      ! thermal velocity for H2SO4 in air (m/s)
      ! https://en.wikipedia.org/wiki/Thermal_velocity
      th = sqrt(8.0_r8*boltz*t0/(pi*molecularWeight *aunit))

      ! calculating microphysical parameters from equations in Ch. 8 of Seinfeld & Pandis (1998):
      ! mean free path for molec in air (m)
      mfv = 1.0_r8/(pi*sqrt(1.0_r8+MolecularWeight/Mair)*(radair+radmol)**2*p0/(boltz*t0))

      ! Solve eqn 11-4.4 in Poling et al
      ! (A bit hard to follow units here, but result in the book is in cm2/s)..
      ! so scale by "cm2Tom2" to get m2/sec
! TODO: what is Vad?
      diff = cm2Tom2   &
         *0.00143_r8*t0**1.75_r8     &
         /((p0/1.0e5_r8)*sqrt(Mdual)   &
         *(((Vad)**aThird+(Vadair)**aThird)**2))

      do ibin = 1, aero_props%nbins()         !all bins receive condensation
          ! Correct for non-continuum effects, formula is from
          ! Chuang and Penner, Tellus, 1995, sticking coeffient from
          ! Vignati et al, JGR, 2004
          diff_coeff(ibin) = diff  &    !original diffusion coefficient
               / ( bin_centers(ibin) / (bin_centers(ibin) + mfv )  &  ! non-continuum correction factor
               +4.0_r8*diff/ (1._r8*th*bin_centers(ibin) ) )
      enddo

      !Find sink per particle in mode "imode"
      !Eqn 13 in Kulmala et al, Tellus 53B, 2001, pp 479
      !http://onlinelibrary.wiley.com/doi/10.1034/j.1600-0889.2001.530411.x/abstract

      !smb++ sectional
      cond_sink_norm = 0.0_r8
      do ibin = 1, aero_props%nbins()
         ! Since we do not sum over bins, it is slightly different than above (normnk=1, no summing)
         cond_sink_norm(ibin) =  &
                                                + 4.0_r8*pi                                    &
                                                * diff_coeff(ibin) &    ![m2/s] diffusion coefficient
                                                * bin_centers(ibin)     ![m] radius of bin
      end do

      !Initialize output
      call phys_getopts(history_aerosol_out = history_aerosol)

  end subroutine condensation_init

  subroutine condtend_sub_super(lchnk,  q, cond_vap_gasprod, temperature, &
               pmid, pdel, dt, ncol, pblh,zm,qh20, aero_props, aero_state)
      ! Calculates nucleation rate and condensation rate of aerosols
      !
      ! This method calls the condtend method. If the timestep needs to be split in two,
      ! this will be done here and the condtend method will be called several times.
      ! This method also writes output once condend is done.
      !

      use cam_history, only: outfld,fieldname_len
      !++smb: add coagulation for npf:
      !use koagsub,    only: normalizedcoagulationsink,receivermode,numberofcoagulationreceivers ! h2so4 and soa nucleation(cka)      !--smb: add coagulation for npf:

      use constituents,    only: pcnst  ! h2so4 and soa nucleation (cka)

      implicit none

      type(sectional_aerosol_properties), intent(in) :: aero_props
      type(sectional_aerosol_state), intent(in) :: aero_state


      ! arguments
      integer,  intent(in)    :: lchnk                      ! chunk identifier
      integer,  intent(in)    :: ncol                       ! number of columns
      real(r8), intent(in)    :: temperature(:,:)    ! Temperature (K)
      real(r8), intent(in)    :: pmid(:,:)           ! [Pa] pressure at mid point
      real(r8), intent(in)    :: pdel(:,:)           ! [Pa] difference in grid cell
      real(r8), intent(in)    :: cond_vap_gasprod(:,:,:) ! TMR [kg/kg/sec]] production rate of H2SO4 (gas prod - aq phase uptake)
      real(r8), intent(in)    :: dt                         ! Time step
      ! Needed for soa nucleation treatment
      real(r8), intent(in)    :: pblh(:)               ! pbl height (m)
      real(r8), intent(in)    :: zm(:,:)           ! midlayer geopotential height above the surface (m) (pver+1)
      real(r8), intent(in)    :: qh20(:,:)          ! specific humidity (kg/kg)

      real(r8), intent(inout) :: q(:,:,:) ! TMR [kg/kg] including moisture


      real(r8) :: q_t0(pcols,pver,gas_pcnst) ! mass before subroutine.

         logical                        :: history_aerosol
         character(128)                 :: long_name
         character(8)                   :: unit

      real(r8)             :: dt_local

      !output:
      real(r8), dimension(pcols, gas_pcnst)            :: coltend
      real(r8), dimension(pcols, gas_pcnst)            :: coltend_dummy
      real(r8) :: nuclrate_pbl(pcols,pver) ![kg/kg] tracer lost
      real(r8) :: nuclrate(pcols,pver) ![kg/kg] tracer lost
      real(r8) :: formrate_pbl(pcols,pver) ![kg/kg] tracer lost
      real(r8) :: formrate(pcols,pver) ![kg/kg] tracer lost
      real(r8) :: h2so4nucl(pcols,pver) ! h2so4 in nucleation code
      real(r8) :: orgnucl(pcols,pver) ! organics in nucleation code
      real(r8) :: grh2so4(pcols,pver) ! growth rate h2so4
      real(r8) :: grsoa(pcols,pver) ! growth rate SOA
      real(r8) :: coagnucl(pcols,pver) ! coagulation in nucleation
      real(r8), allocatable :: leaveSec(:,:,:) ![kg/kg] tracer lost
      real(r8), allocatable :: leaveSec_dummy(:,:,:) ![kg/kg] tracer lost
      logical  :: notDone ! if not done, continues
      logical  :: split_dt ! whether timestep is split or not
      integer  :: nr_dt, cnt,i,j,k  !number of runs, counter, counter, counter
      real(r8), allocatable :: numconc_old(:,:,:) ![#/m3] number concentration before
      real(r8), allocatable :: numconc_new(:,:,:)![#/m3] number concentration new
      real(r8) :: dummy_nc ![#/m3] number concentration
      integer   :: ibin ! index for bin
      !integer   :: tracerIndex
      real(r8)  :: rhoAir
      character(18) :: fieldname_receiver
       real(r8), pointer :: tmp_num(:,:)

      allocate(numconc_old(ncol, pver, aero_props%nbins()))
      allocate(numconc_new(ncol, pver, aero_props%nbins()))

      !initialization
      numconc_old = 0.0_r8
      numconc_new = 0.0_r8

      q_t0(:,:,:)    = q(:,:,:) ! in case timestep needs to be decreased.
      coltend(:,:)   = 0.0_r8
      coltend_dummy(:,:) = 0.0_r8
      nuclrate_pbl(:,:) = 0.0_r8
      nuclrate(:,:)  = 0.0_r8
      formrate_pbl(:,:) = 0.0_r8
      formrate(:,:)  = 0.0_r8
      h2so4nucl(:,:) = 0.0_r8
      orgnucl(:,:)   = 0.0_r8
      grh2so4(:,:)   = 0.0_r8
      grsoa(:,:)     = 0.0_r8
      coagnucl(:,:)  = 0.0_r8
      dt_local       = dt / 2.0_r8 ! always half timestep
      notDone        = .TRUE.
      split_dt       = .FALSE.
      nr_dt          = 2
      cnt            = 1

      ! get bin number concentrations
      tmp_num => null()
      do ibin = 1, aero_props%nbins()
          call aero_state%get_ambient_num(ibin, tmp_num)
          numconc_old(:,:,ibin) = tmp_num(:,:)
      end do

      ! run until no need to split time any longer
      do while (notDone )

           call condtend_sub(lchnk, q, cond_vap_gasprod, temperature, &
                           nuclrate,nuclrate_pbl,formrate, formrate_pbl, coagnucl, &
                           orgnucl, h2so4nucl,grsoa, grh2so4, &
                           coltend_dummy, split_dt, &
                           leaveSec_dummy,&
                           pmid, pdel, dt_local, ncol, pblh, zm, qh20, aero_props)
           coltend= coltend + coltend_dummy*dt_local ! divides by timestep at end
           leaveSec= leaveSec + leaveSec_dummy*dt_local
           if (split_dt) then ! If split timestep: split timestep
               dt_local = dt_local/2.0_r8
               cnt      = 1
               nr_dt    = nr_dt*2
               q(:,:,:) = q_t0(:,:,:)
               coltend(:,:) = 0.0_r8
               coltend_dummy(:,:) = 0.0_r8
               nuclrate_pbl(:,:) = 0.0_r8
               nuclrate(:,:) = 0.0_r8
               formrate_pbl(:,:) = 0.0_r8
               formrate(:,:)  = 0.0_r8
               h2so4nucl(:,:) = 0.0_r8
               orgnucl(:,:)   = 0.0_r8
               grh2so4(:,:)   = 0.0_r8
               grsoa(:,:)     = 0.0_r8
               coagnucl(:,:)  = 0.0_r8
               notDone        = .TRUE.

           else if (nr_dt .eq. cnt) then ! if not, check if count is eq to number of splits
               !if (nr_dt .eq. cnt) then
               notDone=.FALSE.
           else ! if not done and no need to split timestep again, add one to count.
               cnt=cnt+1
           end if
      end do

      ! divide by timestep:
      leaveSec =   leaveSec/dt
      coltend  =   coltend/dt
      nuclrate =   nuclrate/dt
      nuclrate_pbl = nuclrate_pbl/dt
      formrate =   formrate/dt
      formrate_pbl = formrate_pbl/dt
      h2so4nucl =  h2so4nucl/dt
      orgnucl  =   orgnucl/dt
      grh2so4  =   grh2so4/dt
      grsoa    =   grsoa/dt
      coagnucl =  coagnucl/dt

      ! write output
      call outfld('NUCLRATE', nuclrate, pcols   ,lchnk)
      call outfld('NUCLRATE_pbl', nuclrate_pbl, pcols   ,lchnk)
      call outfld('FORMRATE', formrate, pcols   ,lchnk)
      call outfld('FORMRATE_pbl', formrate_pbl, pcols   ,lchnk)
      call outfld('COAGNUCL', coagnucl, pcols   ,lchnk)
      call outfld('GRH2SO4', grh2so4, pcols   ,lchnk)
      call outfld('GRSOA', grsoa, pcols   ,lchnk)
      call outfld('GR', grsoa+grh2so4, pcols   ,lchnk)
      call outfld('ORGNUCL', orgnucl, pcols, lchnk)
      call outfld('H2SO4NUCL', h2so4nucl, pcols, lchnk)
      call outfld('leaveSecH2SO4', leaveSec(:,:,1), pcols,lchnk)
      call outfld('leaveSecSOA', leaveSec(:,:,2), pcols,lchnk)

      call phys_getopts(history_aerosol_out = history_aerosol)

 !     if(history_aerosol)then

!         long_name=trim(solsym(chemistryIndex(l_so4_a1)))//"condTend"
!         call outfld(long_name, coltend(:ncol,chemistryIndex(l_so4_a1)),pcols,lchnk)
!         long_name=trim(solsym(chemistryIndex(l_soa_a1)))//"condTend"
!         call outfld(long_name, coltend(:ncol,chemistryIndex(l_soa_a1)),pcols,lchnk)
!         long_name=trim(solsym(chemistryIndex(l_so4_na)))//"condTend"
!         call outfld(long_name, coltend(:ncol,chemistryIndex(l_so4_na)),pcols,lchnk)
!         long_name=trim(solsym(chemistryIndex(l_soa_na)))//"condTend"
!         call outfld(long_name, coltend(:ncol,chemistryIndex(l_soa_na)),pcols,lchnk)

         !call aerosect_write2file(q,lchnk,ncol,pmid, temperature)

! TODO: make some nice output
 !        do i = 1, aero_props%nbins()
 !          do j=1, aero_props%nspecies_tot()
 !              WRITE(long_name,'(A,I2.2,A)') trim(secSpecNames(j)),i,'_condTend'
 !              call outfld(long_name, coltend(:ncol, chemistryIndex(secConstIndex(j,i))), pcols,lchnk)
 !          end do !j
 !        end do !i
!         end if

       ! extra output:
      numconc_new(:,:,:)=0.0_r8
      do ibin = 1, aero_props%nbins()
          call aero_state%get_ambient_num(ibin, tmp_num)
          numconc_new(:,:,ibin) = tmp_num(:,:)
      end do

      do ibin=1,aero_props%nbins()
           WRITE(fieldname_receiver,'(A,I2.2,A)') 'nrSEC', ibin,'_diff'
           call outfld(trim(fieldname_receiver), (numconc_new(:,:,ibin)-numconc_old(:,:,ibin)),  pcols,lchnk)

      end do

      deallocate(numconc_old, numconc_new)

end subroutine condtend_sub_super


   subroutine condtend_sub(lchnk,  q, cond_vap_gasprod, temperature,            &
           !smb++sectional
                nuclrate,nuclrate_pbl_o, formrate, formrate_pbl_o, coagnucl_o,  &
                orgnucl_o, h2so4nucl_o, grsoa_o, grh2so4_o,                     &
                coltend_o, split_dt,                                            &
                leaveSec,                                                       &
           !smb--sectional
               pmid, pdel, dt, ncol, pblh,zm,qh20,                              &
               aero_props, aero_state)

      ! sub method.
      ! calculate the sulphate nucleation rate, and condensation rate of
      ! aerosols used for parameterising the transfer of externally mixed
      ! aitken mode particles into an internal mixture.
      ! note the parameterisation for conversion of externally mixed particles
      !  used the h2so4 lifetime onto the particles, and not a given
      ! increase in particle radius. will be improved in future versions of the model
      ! added input for h2so4 and soa nucleation: soa_lv_gasprod, soa_sv_gasprod, pblh,zi,qh20 (cka)

      use cam_history,     only: outfld,fieldname_len
      !++smb: add coagulation for npf:
      !use koagsub,         only: normalizedcoagulationsink,receivermode,numberofcoagulationreceivers ! h2so4 and soa nucleation(cka)
      !--smb: add coagulation for npf:
      use constituents,    only: pcnst  ! h2so4 and soa nucleation (cka)

      implicit none

      type(sectional_aerosol_properties), intent(in) :: aero_props
      type(sectional_aerosol_state), intent(in) :: aero_state


       !++smb sectional
      real(r8), intent(inout)  :: nuclrate (pcols, pver)            ! Nucleation rate output
      real(r8), intent(inout)  :: nuclrate_pbl_o (pcols, pver)      ! Nucleation rate pbl output
      real(r8), intent(inout)  :: formrate(pcols, pver)             ! Formation rate output
      real(r8), intent(inout)  :: formrate_pbl_o(pcols,pver)        ! Formation rate pbl output

      real(r8), intent(inout)  :: coagnucl_o(pcols, pver)           ! Coagulation sink for npf output
      real(r8), intent(inout)  :: orgnucl_o(pcols, pver)            ! Organics for nucleation output
      real(r8), intent(inout)  :: h2so4nucl_o(pcols, pver)          ! H2SO4 for nucleation output
      real(r8), intent(inout)  :: grsoa_o(pcols, pver)              ! GR from organics output
      real(r8), intent(inout)  :: grh2so4_o(pcols, pver)            ! GR from H2SO4 output
      real(r8), intent(out)    :: coltend_o(pcols, gas_pcnst)       ! column tendency output
      logical,  intent(out)    :: split_dt                          ! if true, time step needs to be split
      real(r8), allocatable    :: leaveSec(:,:,:)                   ![kg/kg] tracer lost
       !--smb sectional


      ! arguments
      integer,  intent(in) :: lchnk                      ! chunk identifier
      integer,  intent(in) :: ncol                       ! number of columns
      real(r8), intent(in) :: temperature(pcols,pver)    ! Temperature (K)
      real(r8), intent(in) :: pmid(pcols,pver)           ! [Pa] pressure at mid point
      real(r8), intent(in) :: pdel(pcols,pver)           ! [Pa] difference in grid cell
      real(r8), intent(inout) :: q(pcols,pver,gas_pcnst) ! TMR [kg/kg] including moisture
      real(r8), intent(in) :: cond_vap_gasprod(pcols,pver) ! TMR [kg/kg/sec]] production rate of H2SO4 (gas prod - aq phase uptake)
      real(r8), intent(in) :: dt                         ! Time step
      ! Needed for soa nucleation treatment
      real(r8), intent(in)    :: pblh(pcols)               ! pbl height (m)
      real(r8), intent(in)    :: zm(pcols,pverp)           ! midlayer geopotential height above the surface (m) (pver+1)
      real(r8), intent(in)    :: qh20(pcols,pver)          ! specific humidity (kg/kg)

      ! local
      character(len=fieldname_len+3) :: fieldname
      integer :: i,k
      integer :: mode_index_donor            ![idx] index of mode donating mass
      integer :: mode_index_receiver         ![idx] index of mode receiving mass
      integer :: tracerIndex
      integer :: l_donor
      integer :: l_receiver
      integer :: iDonor                                 ![idx] counter for externally mixed modes
      !smb++sectional
      real(r8), allocatable :: condensationsink_sec(:)![1/s] loss rate per mode (mixture)
      !smb--sectional
      !smb++sectional
      real(r8), allocatable :: condensationsinkfraction_sec(:,:,:,:) ! [frc]
      !smb--sectional
      real(r8) :: sumCondensationSink(pcols,pver)       ![1/s] sum of condensation sink
      real(r8) :: totalLoss(pcols,pver,gas_pcnst) ![kg/kg] tracer lost
      !smb++sectional
      real(r8), allocatable :: numberconcentration_sec(:,:,:) ![#/m3] number concentration
      !smb--sectional
      real(r8), dimension(pcols, gas_pcnst)            :: coltend

      real(r8), dimension(pcols)                       :: tracer_coltend


      real(r8)       :: intermediateConcentration(pcols,pver)
      real(r8)       :: rhoAir(pcols,pver)                           ![kg/m3] density of air
      ! Volume of added  material from condensate;  surface area of core particle;
      real(r8)       :: volume_shell, area_core,vol_monolayer
      real(r8)       :: frac_transfer                   ! Fraction of hydrophobic material converted to an internally mixed mode
      logical        :: history_aerosol
      character(128) :: long_name                              ![-] needed for diagnostics

      !cka:+
      ! needed for h2so4 and soa nucleation treatment
       integer  :: modeIndexReceiverCoag              ! Index of modes receiving coagulate
       integer  :: iCoagReceiver                      ! counter for species receiving coagulate
       real(r8) :: coagulationSink(pcols,pver)        ! [1/s] coaglation loss for SO4_n and soa_n
        !nuctst3+
        !   real(r8) :: normCSmode1(pcols,pver)       ! normalized coagulation from self coagulation (simplified)
        !nuctst3-
       real(r8), parameter :: lvocfrac=0.5            ! Fraction of organic oxidation products with low enough
                                                      ! volatility to enter nucleation mode particles (1-24 nm)
       real(r8), pointer :: tmp_num(:,:)
       real(r8) :: soa_lv_forNucleation(pcols,pver)   ! [kg/kg] soa gas available for nucleation
       real(r8) :: gasLost(pcols,pver)                ! [kg/kg] budget terms on H2SO4 (gas)
       real(r8) :: fracNucl(pcols,pver)               ! [frc] fraction of gas nucleated
       real(r8) :: firstOrderLossRateNucl(pcols,pver) ! [1/s] first order loss rate due to nucleation
       real(r8) :: nuclso4(pcols,pver)                ! [kg/kg/s] Nucleated so4 mass tendency from RM's parameterization
       real(r8) :: nuclsoa(pcols,pver)                ! [kg/kg/s] Nucleated soa mass tendency from RM's parameterization
       !smb++ sectional
       real(r8) :: dummy  !
       integer  :: ibin ! indices

       !smb-- sectional

       allocate(condensationsink_sec(aero_props%nbins()))
       allocate(numberconcentration_sec(pcols,pver,aero_props%nbins()))
       !Initialize h2so4 and soa nucl variables
       coagulationSink = 0.0_r8
       condensationsinkfraction_sec = 0.0_r8
       numberconcentration_sec = 0.0_r8
       tmp_num => null()

       do ibin = 1, aero_props%nbins()
        ! No looping through species, mmr is added afterwards
          call aero_state%get_ambient_num(ibin, tmp_num)
          numberconcentration_sec(:,:,ibin) = tmp_num(:,:)
       enddo
       do k=1,pver
           do i=1,ncol

                !smb++ sectional
                ! initialize condensation sink for sectional
                condensationSink_sec = 0.0_r8  !Sink to the coming "receiver" of any vapour
                !smb-- sectional

                !NB: The following is duplicated code, coordinate with koagsub!!
                !Initialize number concentration for this receiver

                !Air density
                rhoAir(i,k) = pmid(i,k)/rair/temperature(i,k)

                !smb++ sectional
                ! Calculate number concentration in each bin :
                numberConcentration_sec(i,k, :) = 0.0_r8

                !Go though all bins receiving condensation

                !smb-- sectional


                !smb++ sectional
                ! condensation sink to sectional bin:
                   do ibin = 1, aero_props%nbins()

                      !This is the loss rate a gas molecule will see due to aerosol surface area
                      condensationSink_sec(ibin)   = cond_sink_norm(ibin)  & ![m3/#/s]  per particle
                                                          * numberConcentration_sec(i,k,ibin)             ![#/m3]
                                                          !==> [1/s]
                   end do !Loop over receivers
                !smb-- sectional


                !Find concentration after condensation of all
                !condenseable vapours
                ! smb++ sectional edited:
                ! condensation sink fraction that goes to each bin in sectional scheme
                ! set to zero first:
                condensationSinkFraction_sec(i,k,:,:)=0.0_r8
                ! smb-- sectional

                    !sum of cond. sink for this vapour [1/s]
                    !smb++ sectional
                    ! Need to add condensation sink to sectional scheme for particles in sectional scheme
                    ! However, not all tracers may contribute:
                    ! assumes same order of gasses in sectional and other (1: H2SO4,2: SOA_LV, 3: SOA_SV)
                        sumCondensationSink(i,k) = sum(condensationSink_sec( :))
                        ! Keeps track of the fraction of the condensate to the sectional bins for each tracer
                        condensationSinkFraction_sec(i,k,:) = condensationSink_sec( :) &
                                /(sumCondensationSink(i,k)+1.e-30_r8)![frc]
                    !smb-- sectional


                !Solve the intermediate (end of timestep) concentration using
                !euler backward solution C_{old} + P *dt - L*C_{new}*dt = C_{new} ==>
                !Cnew -Cold = prod - loss ==>
                intermediateConcentration(i,k) = &
                                     ( q(i,k,l_h2so4_chem) + cond_vap_gasprod(i,k)*dt ) &
                                     / (1.0_r8 + sumCondensationSink(i,k)*dt)



                !Assume only a fraction of ORG_LV left can contribute to nucleation


                modeIndexReceiverCoag = 0
                !Sum coagulation sink for nucleated so4 and soa particles over all receivers of coagulate. Needed for RM's nucleation code
                !OBS - looks like RM's coagulation sink is multiplied by 10^-12??


           end do !index i
       end do !index k

       !Calculate nucleated masses of so4 and soa (nuclso4, nuclsoa)
       !following RM's parameterization (cka)
!       call aeronucl(lchnk,ncol,temperature, pmid, qh20, &
!                   intermediateConcentration(:,:), soa_lv_forNucleation, &
!                   coagulationSink, nuclso4, nuclsoa, zm, pblh, &
                   !smb++ sectional
!                   nuclrate,nuclrate_pbl_o, formrate, formrate_pbl_o, &
!                   orgnucl_o, h2so4nucl_o, grsoa_o, grh2so4_o, dt &
!                   ,bin_centers(1) &
                   !smb-- sectional
!               )

       !smb++ sectional
       coagnucl_o(:,:) = coagnucl_o(:,:) + coagulationSink(:,:)*dt
       !smb-- sectional
       firstOrderLossRateNucl = 0.0_r8
       do k=1,pver
          do i=1,ncol

             !First order loss rate (1/s) for nucleation
              !smb++ added check for 0
              if (intermediateConcentration(i,k) .eq. 0) then
                  firstOrderLossRateNucl(i,k)=0._r8
              else
                  firstOrderLossRateNucl(i,k) = nuclSo4(i,k)/intermediateConcentration(i,k)

              end if
! organics here
              !smb-- added check for 0
             !First order loss rate (1/s) for nucleation

                !Solve implicitly (again)
                !C_new - C_old =  PROD_{gas} - CS*C_new*dt - LR_{nucl}*C_new =>
                intermediateConcentration(i,k) = &
                               ( q(i,k,l_h2so4_chem) + cond_vap_gasprod(i,k)*dt ) &
                               / (1.0_r8 + sumCondensationSink(i,k)*dt + firstOrderLossRateNucl(i,k)*dt)

                !fraction nucleated
                fracNucl(i,k) = firstOrderLossRateNucl(i,k) &
                                     /(firstOrderLossRateNucl(i,k) + sumCondensationSink(i,k))
                !From budget, we get: lost = prod -cnew + cold
                gasLost(i,k) = cond_vap_gasprod(i,k)*dt   & !Produced
                                     + q(i,k,l_h2so4_chem)            & !cold
                                     - intermediateConcentration(i,k)    !cnew


             !Add nuceated mass to so4_na mode
             !smb++ sectional don't add to so4_na directly, must go to sectional scheme (later)
             !q(i,k,chemistryIndex(l_so4_na)) =  q(i,k,chemistryIndex(l_so4_na))       &
             !            + gasLost(i,k,COND_VAP_H2SO4)*fracNucl(i,k,COND_VAP_H2SO4)

             !H2SO4 condensate
             do ibin=1, aero_props%nbins()
                     aero_state%bin_numconc(i, k, ibin) =  &
                !    q(i,k,chemistryIndex(secConstIndex(1,ind_sec))) =                               &
                            !q(i,k,chemistryIndex(secConstIndex(1,ind_sec)))                         &
        ! TODO: check if unit is mass or nr
                            aero_state%bin_numconc(i,k,ibin)  &
                            + gasLost(i,k)*(1.0_r8-fracNucl(i,k))     &
                            *condensationSinkFraction_sec(i,k, ibin)  ! fraction to the particular bin
             end do
             !smb-- sectional
             !H2SO4 condensate
!             q(i,k,chemistryIndex(l_so4_a1)) = q(i,k,chemistryIndex(l_so4_a1))         &
!                            + gasLost(i,k)*(1.0_r8-fracNucl(i,k)) &
                            !smb++ sectional must substract the fraction which goes to the sectional particles:
!                            *(1-sum(condensationSinkFraction_sec(i,k,:)))
                            !smb-- sectional

             !Add nucleated mass to soa_na mode
             !smb++sectional sectional don't add to so4_na directly, must go to sectional scheme (done later)
! organics here

             !condenseable vapours
  !           q(i,k,chemistryIndex(l_h2so4))  = intermediateConcentration(i,k)

             !smb++sectional grow particles in sectional scheme:
    ! TODO: make Sec_movemass routine
    ! removed median radius because not used
             call sec_moveMass(q(i,k,:), numberConcentration_sec(i,k,:), leaveSec(i,k,:), &
                            rhoAir(i,k), split_dt, aero_props)
             ! Add nucleated mass to first bin of sectional scheme:
             q(i,k,chemistryIndex(secConstIndex(1,1))) =  q(i,k,chemistryIndex(secConstIndex(1,1)))       &
                         + gasLost(i,k)*fracNucl(i,k)

            ! Add mass from sectional scheme to so4_na and soa_na:
!             q(i,k,chemistryIndex(l_so4_na)) = q(i,k,chemistryIndex(l_so4_na))         &
!                                    +leaveSec(i,k, 1)

             !smb--sectional

          end do !physical index k
       end do    !physical index i

       !Output for diagnostics
       call phys_getopts(history_aerosol_out = history_aerosol)

       if(history_aerosol)then
          coltend(:ncol,:) = 0.0_r8

          !smb++ sectional
          ! Remove so4_n ---> directly into so4_na
          coltend(:ncol,chemistryIndex(secConstIndex(1,1))) = coltend(:ncol,chemistryIndex(secConstIndex(1,1))) + &
          !smb-- sectional
                                                 sum(                                         &
                                                    gasLost(:ncol,:)           &
                                                    *fracNucl(:ncol,:)*pdel(:ncol,:) , 2 &
                                                    )/gravit/dt

          ! Remove so4_n ---> directly into so4_na
          !smb++ put in how much leaves sectional scheme:
!          coltend(:ncol,chemistryIndex(l_so4_na)) = coltend(:ncol,chemistryIndex(l_so4_na)) + &
!                                                 sum(                                         &
!                                                    leaveSec(:ncol,:)           &
!                                                    *pdel(:ncol,:) , 2 &
!                                                    )/gravit/dt
          !smb--sectional
          !Take into account H2SO4 (gas) condensed in budget
!          coltend(:ncol,chemistryIndex(l_so4_a1)) = coltend(:ncol,chemistryIndex(l_so4_a1)) + &
!                                                 sum(                                         &
!                                                    gasLost(:ncol,:)           &
                                                    !smb++sectional subtract fraction to sectional scheme
!                                                    *(1-sum(condensationSinkFraction_sec(:ncol,:,:),3)) &
                                                    !smb--sectional
!                                                    *(1.0_r8 - fracNucl(:ncol,:))*pdel(:ncol,:) , 2 &
!                                                    )/gravit/dt

          !Take into account soa_lv (gas) nucleated in budget
          !smb++ sectional
! organics here
          !smb++ sectional: putt in how condenses on sectional:
          do ibin=1,aero_props%nbins()
               coltend(:ncol,chemistryIndex(secConstIndex(1,ibin))) = coltend(:ncol, chemistryIndex(secConstIndex(1,ibin)))+ &
                                                    sum(                    &
                                                    gasLost(:ncol, :)    &
                                                    *(condensationSinkFraction_sec(:ncol,:,ibin)) &
                                                    *(1.0_r8 - fracNucl(:ncol, :))*pdel(:ncol,:),2) &
                                                    /gravit/dt
          end do
! organics here

          coltend_o(:,:)=coltend(:,:)

       endif


       return
   end subroutine condtend_sub

subroutine sec_moveMass(massDistrib, numberConc_old, leave_sec, rhoAir, decrease_dt, aero_props)
    ! Moves tracer mass from on bin to the other based on condensational/coagulation growth.
    ! Based on Jacobson Fundamentals of Atmospheric Modeling, second edition (2005),
    ! Chapter   13.5
!    use aerosoldef, only : chemistryIndex
    class(sectional_aerosol_properties), intent(in) :: aero_props

    real(r8), intent(in)    :: numberConc_old(:)    ! numbr concentration before growth
    real(r8), intent(inout) :: massDistrib(:)       ! mass in each tracer
    real(r8), intent(out)   :: leave_sec(:)         ! the mass that leaves sectional scheme
    logical,  intent(out)   :: decrease_dt          ! if set to True, time step is divided
                                                    ! and the procedure is re run

    real(r8), allocatable   :: numberConc_new(:)    ! number concentration after growth  dimension(secNrSpec, secNrBins)
    real(r8), allocatable   :: volume(:)            ! volume of particle in bin dimension(secNrBins)
    real(r8), allocatable   :: volume_new(:)        ! volume after growth dimension(secNrBins)
    real(r8), allocatable   :: volfrac(:)           ! dimension(secNrSpec,secNrBins)
    real(r8)                :: xfrac                ! fraction to stay in bin

    real(r8), parameter     :: pi = 3.141592654_r8
    real(r8)                :: rhoAir               ! Density of air
    integer                 :: ibin

    allocate(numberConc_new(aero_props%nbins()))
    allocate(volume(aero_props%nbins()))
    allocate(volume_new(aero_props%nbins()))
    allocate(volfrac(aero_props%nbins()))

    return
 !   decrease_dt=.FALSE.
    !compute volume in each bin with condensation (by mass) and by
    !numberconcentration
!    do ibin = 1, aero_props%nbins()
!            volume_new(ibin) = 0.0_r8
!            volfrac(ibin) = 0.0_r8
       !     do indSpec = 1, secNrSpec! calculate volume in each bin by mass/density! m3
!                    if (numberConc_old(ibin)<1.e-30_r8) then
!                            volume_new(ibin)=0.0_r8
!                    else
!                        volume_new(ibin) = volume_new(ibin) + massDistrib(chemistryIndex(secConstIndex(ibin)))/&
!                                rhopart_sec * rhoAir/&
!                                (numberConc_old(ibin))
!                    end if

!                    volfrac(ibin)=massDistrib(chemistryIndex(secConstIndex(ibin)))/&
!                            rhopart_sec*rhoAir
                            !kg/kg(air)*[kg(air)/m3(air)][kg/m3]--> m3/m3(air)
        !    end do ! calculate volume in each bin by numberconcentration (volume from before condenstion)
!            volfrac(ibin)=volfrac(ibin)/(sum(volfrac(ibin))+1.E-50_r8)
!            if (volfrac(ibin)<1.e-50_r8) then
!                    volfrac(ibin)=0.0_r8
!            end if
            ! calculate volume in each bin by mass/density! m3
!            volume(ibin) =  pi * secMeanD(ibin)**3/6._r8
            ! calculate volume in each bin by numberconcentration (volume from before condenstion)
!    end do
!    numberConc_new(:) = 0._r8
!    do ibin =  1, aero_props%nbins()-1
            ! fraction to stay in bin
!            xfrac=(volume(ibin+1)-volume_new(ibin)) &
!                            /(volume(ibin+1)-volume(ibin))
!            if (numberConc_old(ibin)<1.e-30) then
!                    xfrac=1.0_r8
!            end if
!            if (xfrac .le. 0._r8) then      ! if the fraction to stay is equal to
                                            ! less than zero, then the
                                            ! aerosols have grown too large
                                            ! for the next bin and we will
                                            ! want to decrease the time step
                                            ! to avoid this.
!                    decrease_dt=.TRUE.
!            end if

!            if (xfrac .le. 0._r8) then
!                    decrease_dt=.TRUE.
!            end if
!            xfrac=max(0._r8, min(1._r8,xfrac))
      !      do indSpec= 1, secNrSpec
!                    numberConc_new(ibin) = numberConc_new(ibin) + &
!                                    xfrac*numberConc_old(ibin) &
!                                    *volfrac(ibin)
!                    numberConc_new(ibin+1) = numberConc_new(ibin+1) + &
!                                    (1-xfrac)*numberConc_old(ibin)   &
!                                    *volfrac(ibin)


!            end do

!    end do

!    xfrac = (max_diameter**3 * pi/6.0_r8 - volume_new(aero_props%nbins())) &
!                            /(max_diameter**3*pi/6.0_r8-volume(aero_props%nbins()))

    ! if less than or 0 % stays in bin, we must decrease timestep
!    if (xfrac .le. 0._r8) then
!            decrease_dt=.TRUE.
!    end if

!    xfrac=max(0._r8, min(1._r8,xfrac))

   ! do indSpec=1, secNrSpec
!            numberConc_new(aero_props%nbins()) = numberConc_new(aero_props%nbins()) + &
!                                            xfrac * numberConc_old(aero_props%nbins()) &
!                                            * volfrac(aero_props%nbins())
!            leave_sec = & !massDistrib(chemistryIndex(secConstIndex(indSpec, aero_props%nbins())))*(1-xfrac)
!                    pi * max_diameter**3 / 6.0_r8 * rhopart_sec/rhoAir &   ! [m3_aer/#]*[kg_aer/m3_aer]/[kg_air/m3_air]--> [kg_aer/kg_air/#][m3_air]
!                                            * (1-xfrac) * numberConc_old(aero_props%nbins()) &          ! *[#/m3_air] --> kg_aer/kg_air
!                                            * volfrac(aero_props%nbins())
   ! end do
!    do ibin=1,aero_props%nbins()
           ! do indSpec=1,secNrSpec !Assume
!                    massDistrib(chemistryIndex(secConstIndex,ibin)) = &! &!massDistrib(secConstIndex(indSpec,ibin))+&
!                            rhopart_sec/rhoAir &!* massfrac(indSpec,ibin)* numberConc_new(ibin)! &
!                            * numberConc_new(ibin) * pi * secMeanD(ibin)**3/6.0_r8 !&
           ! end do
!    end do



end subroutine sec_moveMass


end module condtend