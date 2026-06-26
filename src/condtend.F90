module condtend

   use phys_control, only: phys_getopts
   use chem_mods,    only: gas_pcnst
   use mo_tracname,  only: solsym
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use const
   use cam_history,  only: outfld
   use aerosoldef
   use physconst,    only: rair, gravit, pi
   use commondefinitions
   use chem_mods, only: adv_mass !molecular weights from mozart
   !smb++ sectional
   use aero_sectional,    only: secNrSpec,secNrBins, secMeanD !
   !smb-- sectional
!soa

   save

   integer, parameter :: N_COND_VAP = 3
   integer, parameter :: COND_VAP_H2SO4 = 1
   integer, parameter :: COND_VAP_ORG_LV = 2
   integer, parameter :: COND_VAP_ORG_SV = 3

   real(r8), public, dimension(0:nmodes,N_COND_VAP)  :: normalizedCondensationSink       ![m3/#/s] condensation sink per particle in mode i
   ! smb++ sectional
   real(r8), public, dimension(secNrSpec, secNrBins)  :: normalizedCondensationSink_sec       ![m3/#/s] condensation sink per particle in bin i
   ! smb-- sectional

   integer, private, dimension(gas_pcnst) :: lifeCycleReceiver                ! [-] array of transformation of life cycle tracers
   real(r8), private, dimension(0:nmodes,N_COND_VAP) :: stickingCoefficient              ! [-] stickingCoefficient for H2SO4 on a mode
   integer, private, dimension(N_COND_VAP) :: cond_vap_map

! Assumed number of monolayers
  real(r8), parameter, private :: n_so4_monolayers_age = 3.0_r8

  real(r8), parameter, public :: &
              dr_so4_monolayers_age = n_so4_monolayers_age * 4.76e-10_r8
! thickness of the so4 monolayers (m)
! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3 as in MAM


contains

   subroutine registerCondensation()

      implicit none


      integer                        :: iDonor
      integer                        :: l_donor
      integer                        :: tracerIndex
      integer                        :: mode_index_donor

      !These are the lifecycle-species which receive mass when
      !the externally mixed modes receive condensate,
      !e.g. the receiver of l_so4_n mass is the tracer l_so4_na
      lifeCycleReceiver(:) = -99
      lifeCycleReceiver(chemistryIndex(l_bc_n))   = chemistryIndex(l_bc_a)    !create bc int mix from bc in mode 12
      lifeCycleReceiver(chemistryIndex(l_bc_ni))  = chemistryIndex(l_bc_ai)   !create bc int mix from bc in mode 14
      lifeCycleReceiver(chemistryIndex(l_om_ni))  = chemistryIndex(l_om_ai)
      !!create om int mix from om in mode 14
      lifeCycleReceiver(chemistryIndex(l_bc_ax))  = chemistryIndex(l_bc_ai)
      !!create bc int mix from bc in mode 0. Note Mass is conserved but not number

      !Sticking coeffcients for H2SO4 condensation
      !See table 1 in Kirkevag et al (2013)
      !http://www.geosci-model-dev.net/6/207/2013/gmd-6-207-2013.html
      !Note: In NorESM1, sticking coefficients of the externally mixed modes were
      !used for the internally mixed modes in modallapp. In condtend the internally
      !mixed modes had sticking coefficient = 1.0
      !This might be correct, but is too confusing, so here just
      !assign based on background aerosol and table 1 in Kirkevag et al
      stickingCoefficient(:,:) = 1.0_r8
      stickingCoefficient(MODE_IDX_BC_EXT_AC,:) = 0.3_r8
      stickingCoefficient(MODE_IDX_BC_AIT,:) = 0.3_r8
      stickingCoefficient(MODE_IDX_OMBC_INTMIX_COAT_AIT,:) = 0.5_r8
      stickingCoefficient(MODE_IDX_DST_A2,:) = 0.3_r8
      stickingCoefficient(MODE_IDX_DST_A3,:) = 0.3_r8
      stickingCoefficient(MODE_IDX_BC_NUC,:) = 0.3_r8
      stickingCoefficient(MODE_IDX_OMBC_INTMIX_AIT,:) = 0.5_r8


   end subroutine registerCondensation

!===============================================================================

   subroutine initializeCondensation()

      !condensation coefficients:
      !Theory: Poling et al, "The properties of gases and liquids"
      !5th edition, eqn 11-4-4

      use cam_history,     only: addfld, add_default, fieldname_len, horiz_only
      !smb++ sectional
      use aero_sectional,     only: secSpecNames, secNrSpec, secNrBins, secMeanD !smb:sectional
      !smb-- sectional

      implicit none

      real(r8), parameter :: aunit = 1.6606e-27_r8  ![kg] Atomic mass unit
      real(r8), parameter :: boltz = 1.3806e-23_r8   ![J/K/molec]
      real(r8), parameter :: t0 = 273.15_r8         ![K] standard temperature
      real(r8), parameter :: p0 = 101325.0_r8       ! [Pa] Standard pressure
      real(r8), parameter :: radair = 1.73e-10_r8   ![m] Typical air molecule collision radius
      real(r8), parameter :: Mair = 28.97_r8        ![amu/molec] Molecular weight for dry air
      !Diffusion volumes for simple molecules [Poling et al], table 11-1
      real(r8), dimension(N_COND_VAP), parameter :: vad = (/51.96_r8, 208.18_r8, 208.18_r8/) ![cm3/mol]
      real(r8), parameter :: vadAir       = 19.7_r8                                          ![cm3/mol]
      real(r8), parameter :: aThird = 1.0_r8/3.0_r8
      real(r8), parameter :: cm2Tom2 = 1.e-4_r8       !convert from cm2 ==> m2

      real(r8), dimension(0:100,0:nmodes,N_COND_VAP) :: DiffusionCoefficient   ! [m2/s] Diffusion coefficient
      !smb++ sectional
      real(r8), dimension(secNrSpec, secNrBins) :: DiffusionCoefficientSec   ! [m2/s] Diffusion coefficient sectional
      !smb-- sectional
      character(len=fieldname_len+3) :: fieldname_donor
      character(len=fieldname_len+3) :: fieldname_receiver
      character(128)                 :: long_name
      character(8)                   :: unit

      integer                        :: nsiz !counter for aerotab sizes
      integer                        :: iChem             !counter for chemical species
      integer                        :: mode_index_donor  !index for mode
      integer                        :: iMode             !Counter for mode
      !smb++ sectional
      integer                        :: iBin             !Counter for bin
      !smb-- sectional
      integer                        :: tracerIndex       !counter for chem. spec

      logical                        :: history_aerosol
      logical                        :: isAlreadyOnList(gas_pcnst)
      integer                        :: cond_vap_idx

      real(r8), dimension(N_COND_VAP) :: mfv  ![m] mean free path
      real(r8), dimension(N_COND_VAP) :: diff ![m2/s] diffusion coefficient for cond. vap
      real(r8) :: molecularWeight !amu/molec molecular weight
      real(r8) :: Mdual ![molec/amu] 1/M_1 + 1/M_2
      real(r8) :: rho   ![kg/m3] density of component in question
      real(r8) :: radmol ![m] radius molecule
      real(r8), dimension(N_COND_VAP) :: th     !thermal velocity
      !smb++ sectional
      integer                         :: indSpec, indBin ! indexes
      !smb-- sectional

      !Couple the condenseable vapours to chemical species for properties and indexes
      cond_vap_map(COND_VAP_H2SO4) = chemistryIndex(l_h2so4)
      cond_vap_map(COND_VAP_ORG_LV) = chemistryIndex(l_soa_lv)
      cond_vap_map(COND_VAP_ORG_SV) = chemistryIndex(l_soa_sv)

      do cond_vap_idx = 1, N_COND_VAP

         rho = rhopart(physicsIndex(cond_vap_map(cond_vap_idx))) !pick up densities from aerosoldef

         molecularWeight=adv_mass(cond_vap_map(cond_vap_idx))    !pick up molecular weights from mozart

         !https://en.wikipedia.org/wiki/Thermal_velocity
         th(cond_vap_idx) = sqrt(8.0_r8*boltz*t0/(pi*molecularweight*aunit))   ! thermal velocity for H2SO4 in air (m/s)

         !Radius of molecul (straight forward assuming spherical)
         radmol=(3.0_r8*molecularWeight*aunit/(4.0_r8*pi*rho))**aThird    ! molecule radius

         Mdual=2.0_r8/(1.0_r8/Mair+1.0_r8/molecularWeight) !factor of [1/m_1 + 1_m2]

         !calculating microphysical parameters from equations in Ch. 8 of Seinfeld & Pandis (1998):
         mfv(cond_vap_idx)=1.0_r8/(pi*sqrt(1.0_r8+MolecularWeight/Mair)*(radair+radmol)**2*p0/(boltz*t0)) ! mean free path for molec in air (m)

         !Solve eqn 11-4.4 in Poling et al
         !(A bit hard to follow units here, but result in the book is in cm2/s)..
         !so scale by "cm2Tom2" to get m2/sec
         diff(cond_vap_idx) = cm2Tom2   &
            *0.00143_r8*t0**1.75_r8     &
          /((p0/1.0e5_r8)*sqrt(Mdual)   &
          *(((Vad(cond_vap_idx))**aThird+(Vadair)**aThird)**2))

         !Values used in noresm1:
         !real(r8), parameter :: diff = 9.5e-6    !m2/s  diffusion coefficient (H2SO4)
         !real(r8), parameter :: th   = 243.0_r8  !m/s   thermal velocity (H2SO4)
         !real(r8), parameter :: mfv  = 1.65e-8   !m     mean free path (H2SO4)

         !Check values obtained here (H2SO4 / SOA)
         !write(*,*) 'mfv =   ', mfv(cond_vap_idx)     !2.800830854409093E-008 / 1.633546464678737E-008
         !write(*,*) ' diff = ', diff(cond_vap_idx)    !->  9.360361706957621E-006 / !->  4.185923463242946E-006
         !write(*,*) ' th = ', th                      !->  242.818542922924 / 185.421069430852
      end do

      do cond_vap_idx = 1, N_COND_VAP
         do imode = 0, nmodes         !all modes receive condensation
            do nsiz = 1, nBinsTab     !aerotab sizes
            !Correct for non-continuum effects, formula is from
            !Chuang and Penner, Tellus, 1995, sticking coeffient from
            !Vignati et al, JGR, 2004
            !fxm: make "diff ==> diff (cond_vap_idx)
            DiffusionCoefficient(nsiz,imode,cond_vap_idx) = diff(cond_vap_idx)  &    !original diffusion coefficient
               /(                                    &
                  rBinMidPoint(nsiz)/(rBinMidPoint(nsiz)+mfv(cond_vap_idx))  &  !non-continuum correction factor
                  +4.0_r8*diff(cond_vap_idx)/(stickingCoefficient(imode,cond_vap_idx)*th(cond_vap_idx)*rBinMidPoint(nsiz)) &
                 )
            enddo
         end do !receiver modes
      end do
      !smb++ sectional
      do indSpec = 1,secNrSpec ! loop through species
         do indBin = 1, secNrBins         !all bins receive condensation
             !Correct for non-continuum effects, formula is from
             !Chuang and Penner, Tellus, 1995, sticking coeffient from
             !Vignati et al, JGR, 2004
             !fxm: make "diff ==> diff (cond_vap_idx)
             DiffusionCoefficientSec(indSpec, indBin) = diff( indSpec)  &    !original diffusion coefficient
               /(                                    &
                  secMeanD(indBin)*0.5_r8/(secMeanD(indBin)*0.5_r8 + mfv(indSpec))  &  !non-continuum correction factor
                  +4.0_r8*diff(indSpec)/(1._r8*th(indSpec)*0.5_r8*secMeanD(indBin)) &  ! sticking coefficient set to 1
                 )
         enddo
       enddo
      !smb-- sectional

      normalizedCondensationSink(:,:) = 0.0_r8
      !Find sink per particle in mode "imode"
      !Eqn 13 in Kulmala et al, Tellus 53B, 2001, pp 479
      !http://onlinelibrary.wiley.com/doi/10.1034/j.1600-0889.2001.530411.x/abstract
      do cond_vap_idx =1, N_COND_VAP
         do imode = 0, nmodes
            do nsiz = 1, nBinsTab
               normalizedCondensationSink(imode,cond_vap_idx) =  &
                                                normalizedCondensationSink(imode,cond_vap_idx)  &
                                                + 4.0_r8*pi                                    &
                                                * DiffusionCoefficient(nsiz,imode,cond_vap_idx) &    ![m2/s] diffusion coefficient
                                                * rBinMidPoint(nsiz)                &    ![m] look up table radius
                                                * normnk(imode,nsiz)                     ![frc]
            end do
         end do
      end do
      !smb++ sectional
      normalizedCondensationSink_sec(:,:) =0.0_r8
      do indSpec =1, secNrSpec
         do indBin = 1, secNrBins
             ! Since we do not sum over bins, it is slightly different than above (normnk=1, no summing)
               normalizedCondensationSink_sec(indSpec, indBin) =  &
                                                + 4.0_r8*pi                                    &
                                                * DiffusionCoefficientSec( indSpec, indBin) &    ![m2/s] diffusion coefficient
                                                * secMeanD(indBin)*0.5_r8                    ![m] radius of bin
        end do
      end do

      !smb-- sectional

      !Initialize output
      call phys_getopts(history_aerosol_out = history_aerosol)

      isAlreadyOnList(:) = .FALSE.
      do iChem = 1,gas_pcnst
         !Does this tracer have a receiver? If yes: It participate in condensation tendencies
         if(lifeCycleReceiver(iChem) .gt. 0)then
            unit = "kg/m2/s"
            fieldname_donor = trim(solsym(iChem))//"condTend"
            fieldname_receiver = trim(solsym(lifeCycleReceiver(iChem)))//"condTend"
            if(.not. isAlreadyOnList(lifeCycleReceiver(iChem)))then
               call addfld( fieldname_receiver, horiz_only, "A", unit, "condensation tendency" )
               isAlreadyOnList(lifeCycleReceiver(iChem))=.TRUE.
            end if
            call addfld( fieldname_donor, horiz_only, "A", unit, "condensation tendency" )
            if(history_aerosol)then
               call add_default( fieldname_receiver, 1, ' ' )
               call add_default( fieldname_donor   , 1, ' ')
            end if
         end if
      end do
      !Need to add so4_a1, soa_na, so4_na, soa_a1 also (which are not parts of the donor-receiver stuff)
      fieldname_receiver = trim(solsym(chemistryIndex(l_so4_a1)))//"condTend"
      call addfld( fieldname_receiver, horiz_only, 'A', unit, "condensation tendency")
      if(history_aerosol)then
         call add_default( fieldname_receiver, 1, ' ' )
      end if
      fieldname_receiver = trim(solsym(chemistryIndex(l_soa_a1)))//"condTend"
      call addfld( fieldname_receiver, horiz_only, "A", unit, "condensation tendency" )
      if(history_aerosol)then
         call add_default( fieldname_receiver, 1, ' ' )
      end if
      fieldname_receiver = trim(solsym(chemistryIndex(l_so4_na)))//"condTend"
      call addfld( fieldname_receiver, horiz_only, 'A', unit , "condensation tendency" )
      if(history_aerosol)then
         call add_default( fieldname_receiver, 1, ' ' )
      end if
      fieldname_receiver = trim(solsym(chemistryIndex(l_soa_na)))//"condTend"
      call addfld( fieldname_receiver, horiz_only, 'A', unit, "condensation tendency" )
      if(history_aerosol)then
         call add_default( fieldname_receiver, 1, ' ' )
      end if

  !++smb sectional:
   do iChem = 1, secNrSpec
      do imode = 1, secNrBins
        WRITE(fieldname_receiver,'(A,I2.2,A)') trim(secSpecNames(iChem)),imode,'_condTend'
        call addfld(fieldname_receiver, horiz_only, "A", unit,'condensation tendency')
        if (history_aerosol) then
          call add_default(fieldname_receiver, 1, ' ' )
        end if
      end do !imod
   end do !iChem
       ! +++smb extra output smb #todo remove this.
    do imode=1,secNrBins
       WRITE(fieldname_receiver,'(A,I2.2,A)') 'nrSEC', imode,'_diff'
       call addfld(trim(fieldname_receiver),  (/'lev'/), 'A','nr/m3','difference in sectional scheme before and after condensation')
    end do
    !smb-- sectional



   end subroutine initializeCondensation

  subroutine condtend(lchnk,  q, cond_vap_gasprod, temperature, &
               pmid, pdel, dt, ncol, pblh,zm,qh20)
      ! Calculates nucleation rate and condensation rate of aerosols
      !
      ! This method calls the condtend method. If the timestep needs to be split in two,
      ! this will be done here and the condtend method will be called several times.
      ! This method also writes output once condend is done.
      !

      use cam_history,     only: outfld,fieldname_len
      !++smb: add coagulation for npf:
      !use koagsub,         only: normalizedcoagulationsink,receivermode,numberofcoagulationreceivers ! h2so4 and soa nucleation(cka)
      use koagsub,         only: normalizedcoagulationsinknpf,receivermodenpf,numberofcoagulationreceiversnpf ! h2so4 and soa nucleation(cka)
      !--smb: add coagulation for npf:

      use aero_sectional,     only: secSpecNames, secConstIndex, secNrBins, secNrSpec, sec_numberconc

      use constituents,    only: pcnst  ! h2so4 and soa nucleation (cka)

      implicit none

      ! arguments
      integer,  intent(in) :: lchnk                      ! chunk identifier
      integer,  intent(in) :: ncol                       ! number of columns
      real(r8), intent(in) :: temperature(pcols,pver)    ! Temperature (K)
      real(r8), intent(in) :: pmid(pcols,pver)           ! [Pa] pressure at mid point
      real(r8), intent(in) :: pdel(pcols,pver)           ! [Pa] difference in grid cell
      real(r8), intent(inout) :: q(pcols,pver,gas_pcnst) ! TMR [kg/kg] including moisture
      real(r8), intent(in) :: cond_vap_gasprod(pcols,pver,N_COND_VAP) ! TMR [kg/kg/sec]] production rate of H2SO4 (gas prod - aq phase uptake)
      real(r8), intent(in) :: dt                         ! Time step
      ! Needed for soa nucleation treatment
      real(r8), intent(in)    :: pblh(pcols)               ! pbl height (m)
      real(r8), intent(in)    :: zm(pcols,pverp)           ! midlayer geopotential height above the surface (m) (pver+1)
      real(r8), intent(in)    :: qh20(pcols,pver)          ! specific humidity (kg/kg)

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
      real(r8) :: leaveSec(pcols,pver,secNrSpec) ![kg/kg] tracer lost
      real(r8) :: leaveSec_dummy(pcols,pver, secNrSpec) ![kg/kg] tracer lost
      logical  :: notDone ! if not done, continues
      logical  :: split_dt ! whether timestep is split or not
      integer  :: nr_dt, cnt,i,j,k  !number of runs, counter, counter, counter
      real(r8) :: numberconcentration_sec_old(pcols,pver, secnrbins) ![#/m3] number concentration before
      real(r8) :: numberconcentration_sec_new(pcols,pver, secnrbins) ![#/m3] number concentration new
      real(r8) :: dummy_nc ![#/m3] number concentration
      integer   :: indBin ! index for bin
      integer   :: ind_sec ! chemistry index sectional
      !integer   :: tracerIndex
      integer   :: indSpec ! index for condensing volatile species
      real(r8)  :: rhoAir
      character(18) :: fieldname_receiver
      ! LOOP OVER i, k!
      numberconcentration_sec_old(:,:,:)=0.0_r8
      do k=1,pver
         do i=1,ncol
             rhoAir = pmid(i,k)/rair/temperature(i,k)
             do indBin = 1, secNrBins
               !Go through all species in bin
               do indSpec = 1,secNrSpec
                   ind_sec=chemistryIndex(secConstIndex(indSpec, indBin))
                   call sec_numberConc(q(i,k,ind_sec),indSpec, &!SpecNames(tracerIndex), &
                           indBin,rhoAir, &
                           dummy_nc)
                   ! add concentration from species
                   numberConcentration_sec_old(i,k, indBin)=numberconcentration_sec_old(i,k, indBin)+&
                                                                   dummy_nc
               end do
           end do
         end do
      end do


      !initialization
      q_t0(:,:,:)=q(:,:,:) ! in case timestep needs to be decreased.
      coltend(:,:)=0.0_r8
      coltend_dummy(:,:)=0.0_r8
      nuclrate_pbl(:,:)=0.0_r8
      nuclrate(:,:)=0.0_r8
      formrate_pbl(:,:)=0.0_r8
      formrate(:,:)=0.0_r8
      h2so4nucl(:,:)=0.0_r8
      orgnucl(:,:)=0.0_r8
      grh2so4(:,:)=0.0_r8
      grsoa(:,:)=0.0_r8
      coagnucl(:,:)=0.0_r8
      dt_local=dt / 2.0_r8 ! always half timestep
      notDone=.TRUE.
      split_dt=.FALSE.
      nr_dt=2
      cnt=1
      ! run until no need to split time any longer
      do while (notDone )

           call condtend_sub(lchnk, q, cond_vap_gasprod, temperature, &
                           nuclrate,nuclrate_pbl,formrate, formrate_pbl, coagnucl, &
                           orgnucl, h2so4nucl,grsoa, grh2so4, &
                           coltend_dummy, split_dt, &
                           leaveSec_dummy,&
                           pmid, pdel, dt_local, ncol, pblh, zm, qh20)
           coltend= coltend + coltend_dummy*dt_local ! divides by timestep at end
           leaveSec= leaveSec + leaveSec_dummy*dt_local
           if (split_dt) then ! If split timestep: split timestep
               dt_local=dt_local/2.0_r8
               cnt=1
               nr_dt=nr_dt*2
               q(:,:,:)=q_t0(:,:,:)
               coltend(:,:)=0.0_r8
               coltend_dummy(:,:)=0.0_r8
               nuclrate_pbl(:,:)=0.0_r8
               nuclrate(:,:)=0.0_r8
               formrate_pbl(:,:)=0.0_r8
               formrate(:,:)=0.0_r8
               h2so4nucl(:,:)=0.0_r8
               orgnucl(:,:)=0.0_r8
               grh2so4(:,:)=0.0_r8
               grsoa(:,:)=0.0_r8
               coagnucl(:,:)=0.0_r8
               notDone=.TRUE.

           else if (nr_dt .eq. cnt) then ! if not, check if count is eq to number of splits
               !if (nr_dt .eq. cnt) then
               notDone=.FALSE.
           else ! if not done and no need to split timestep again, add one to count.
               cnt=cnt+1
           end if
      end do

      ! divide by timestep:
      leaveSec= leaveSec/dt
      coltend=coltend/dt
      nuclrate=nuclrate/dt
      nuclrate_pbl=nuclrate_pbl/dt
      formrate=formrate/dt
      formrate_pbl=formrate_pbl/dt
      h2so4nucl=h2so4nucl/dt
      orgnucl=orgnucl/dt
      grh2so4=grh2so4/dt
      grsoa=grsoa/dt
      coagnucl=coagnucl/dt
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

      if(history_aerosol)then

         do i=1,gas_pcnst
            if(lifeCycleReceiver(i) .gt. 0 )then
               long_name= trim(solsym(i))//"condTend"
               call outfld(long_name, coltend(:ncol,i), pcols, lchnk)
               long_name= trim(solsym(lifeCycleReceiver(i)))//"condTend"
               call outfld(long_name, coltend(:ncol,lifeCycleReceiver(i)),pcols,lchnk)
            end if
         end do
         long_name=trim(solsym(chemistryIndex(l_so4_a1)))//"condTend"
         call outfld(long_name, coltend(:ncol,chemistryIndex(l_so4_a1)),pcols,lchnk)
         long_name=trim(solsym(chemistryIndex(l_soa_a1)))//"condTend"
         call outfld(long_name, coltend(:ncol,chemistryIndex(l_soa_a1)),pcols,lchnk)
         long_name=trim(solsym(chemistryIndex(l_so4_na)))//"condTend"
         call outfld(long_name, coltend(:ncol,chemistryIndex(l_so4_na)),pcols,lchnk)
         long_name=trim(solsym(chemistryIndex(l_soa_na)))//"condTend"
         call outfld(long_name, coltend(:ncol,chemistryIndex(l_soa_na)),pcols,lchnk)

         !call aerosect_write2file(q,lchnk,ncol,pmid, temperature)
         do i = 1, secNrBins
           do j=1, secNrSpec
               WRITE(long_name,'(A,I2.2,A)') trim(secSpecNames(j)),i,'_condTend'
               call outfld(long_name, coltend(:ncol, chemistryIndex(secConstIndex(j,i))), pcols,lchnk)
           end do !j
         end do !i
         end if

       ! extra output:
      numberConcentration_sec_new(:,:,:)=0.0_r8
      do k=1,pver
         do i=1,ncol
           do indBin = 1, secNrBins
               !Go through all core species in that mode
               do indSpec = 1,secNrSpec
                   ind_sec=chemistryIndex(secConstIndex(indSpec, indBin))
                   rhoAir = pmid(i,k)/rair/temperature(i,k)
                   call sec_numberConc(q(i,k,ind_sec),indSpec, &!SpecNames(indSpec), &
                           indBin,rhoAir, &
                           dummy_nc)
                   numberconcentration_sec_new(i,k, indBin)=numberConcentration_sec_new(i,k, indBin)+&
                                                                   dummy_nc
               end do
           end do
         end do
      end do
      do indBin=1,secNrBins
           WRITE(fieldname_receiver,'(A,I2.2,A)') 'nrSEC', indBin,'_diff'
           call outfld(trim(fieldname_receiver), (numberconcentration_sec_new(:,:,indBin)-numberconcentration_sec_old(:,:,indBin)),  pcols,lchnk)

      end do

end subroutine condtend


   subroutine condtend_sub(lchnk,  q, cond_vap_gasprod, temperature,            &
           !smb++sectional
                nuclrate,nuclrate_pbl_o, formrate, formrate_pbl_o, coagnucl_o,  &
                orgnucl_o, h2so4nucl_o, grsoa_o, grh2so4_o,                     &
                coltend_o, split_dt,                                            &
                leaveSec,                                                       &
           !smb--sectional
               pmid, pdel, dt, ncol, pblh,zm,qh20)

      use aero_sectional,     only: sec_movemass,sec_numberconc, secconstindex,SpecNames, secNrSpec, secnrbins, secmeand !smb:sectional
      use commondefinitions, only: originalnumbermedianradius

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
      use koagsub,         only: normalizedcoagulationsinknpf,receivermodenpf,numberofcoagulationreceiversnpf ! h2so4 and soa nucleation(cka)
      !--smb: add coagulation for npf:
      use constituents,    only: pcnst  ! h2so4 and soa nucleation (cka)

      implicit none
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
      real(r8), dimension(pcols, gas_pcnst) :: coltend_o            ! column tendency output
      logical,  intent(out)    :: split_dt                          ! if true, time step needs to be split
      real(r8), intent(out)    :: leaveSec(pcols, pver, secNrSpec)  ! leaving the sectional scheme
       !--smb sectional


      ! arguments
      integer,  intent(in) :: lchnk                      ! chunk identifier
      integer,  intent(in) :: ncol                       ! number of columns
      real(r8), intent(in) :: temperature(pcols,pver)    ! Temperature (K)
      real(r8), intent(in) :: pmid(pcols,pver)           ! [Pa] pressure at mid point
      real(r8), intent(in) :: pdel(pcols,pver)           ! [Pa] difference in grid cell
      real(r8), intent(inout) :: q(pcols,pver,gas_pcnst) ! TMR [kg/kg] including moisture
      real(r8), intent(in) :: cond_vap_gasprod(pcols,pver,N_COND_VAP) ! TMR [kg/kg/sec]] production rate of H2SO4 (gas prod - aq phase uptake)
      real(r8), intent(in) :: dt                         ! Time step
      ! Needed for soa nucleation treatment
      real(r8), intent(in)    :: pblh(pcols)               ! pbl height (m)
      real(r8), intent(in)    :: zm(pcols,pverp)           ! midlayer geopotential height above the surface (m) (pver+1)
      real(r8), intent(in)    :: qh20(pcols,pver)          ! specific humidity (kg/kg)

      ! local
      character(len=fieldname_len+3) :: fieldname
      integer :: i,k,nsiz
      integer :: mode_index_donor            ![idx] index of mode donating mass
      integer :: mode_index_receiver         ![idx] index of mode receiving mass
      integer :: tracerIndex
      integer :: l_donor
      integer :: l_receiver
      integer :: iDonor                                 ![idx] counter for externally mixed modes
      real(r8) :: condensationSink(0:nmodes, N_COND_VAP)![1/s] loss rate per mode (mixture)
      !smb++sectional
      real(r8) :: condensationsink_sec(secNrSpec, secnrbins )![1/s] loss rate per mode (mixture)
      !smb--sectional
      real(r8) :: condensationSinkFraction(pcols,pver,numberOfExternallyMixedModes,N_COND_VAP) ![frc]
      !smb++sectional
      real(r8) :: condensationsinkfraction_sec(pcols,pver,N_COND_VAP, secnrbins) ![frc]
      !smb--sectional
      real(r8) :: sumCondensationSink(pcols,pver, N_COND_VAP)       ![1/s] sum of condensation sink
      real(r8) :: totalLoss(pcols,pver,gas_pcnst) ![kg/kg] tracer lost
      real(r8) :: numberConcentration(0:nmodes) ![#/m3] number concentration
      !smb++sectional
      real(r8) :: numberconcentration_sec(pcols,pver, secnrbins) ![#/m3] number concentration
      !smb--sectional
      real(r8) :: numberConcentrationExtMix(pcols,pver,numberOfExternallyMixedModes)
      real(r8), dimension(pcols, gas_pcnst)            :: coltend

      real(r8), dimension(pcols)                       :: tracer_coltend


      real(r8) :: intermediateConcentration(pcols,pver,N_COND_VAP)
      real(r8) :: rhoAir(pcols,pver)                           ![kg/m3] density of air
      ! Volume of added  material from condensate;  surface area of core particle;
      real(r8) :: volume_shell, area_core,vol_monolayer
      real (r8) :: frac_transfer                   ! Fraction of hydrophobic material converted to an internally mixed mode
      logical  :: history_aerosol
      character(128)                 :: long_name                              ![-] needed for diagnostics

      !cka:+
      ! needed for h2so4 and soa nucleation treatment
       integer  :: modeIndexReceiverCoag             !Index of modes receiving coagulate
       integer  :: iCoagReceiver                     !counter for species receiving coagulate
       real(r8) :: coagulationSink(pcols,pver)       ![1/s] coaglation loss for SO4_n and soa_n
        !nuctst3+
        !   real(r8) :: normCSmode1(pcols,pver)           !normalized coagulation from self coagulation (simplified)
        !nuctst3-
       real(r8), parameter :: lvocfrac=0.5           !Fraction of organic oxidation products with low enough
                                                      !volatility to enter nucleation mode particles (1-24 nm)
       real(r8) :: soa_lv_forNucleation(pcols,pver)  ![kg/kg] soa gas available for nucleation
       real(r8) :: gasLost(pcols,pver,N_COND_VAP)          ![kg/kg] budget terms on H2SO4 (gas)
       real(r8) :: fracNucl(pcols,pver,N_COND_VAP)               ! [frc] fraction of gas nucleated
       real(r8) :: firstOrderLossRateNucl(pcols,pver,N_COND_VAP) ![1/s] first order loss rate due to nucleation
       real(r8) :: nuclso4(pcols,pver)               ![kg/kg/s] Nucleated so4 mass tendency from RM's parameterization
       real(r8) :: nuclsoa(pcols,pver)               ![kg/kg/s] Nucleated soa mass tendency from RM's parameterization
       integer  :: cond_vap_idx
       !smb++ sectional
       real(r8) :: dummy  !
       integer  :: ind_sec  ! chemistry index sectional tracers
       integer  :: indbin, indSpec ! indices
       !smb-- sectional

       !Initialize h2so4 and soa nucl variables
       coagulationSink(:,:)=0.0_r8
       condensationSinkFraction(:,:,:,:) = 0.0_r8  !Sink to the coming "receiver" of any vapour
       numberConcentrationExtMix(:,:,:) = 0.0_r8
       !smb++ sectional
       condensationsinkfraction_sec(:,:,:, :) = 0.0_r8
       numberconcentration_sec(:,:, :) = 0.0_r8
       !smb-- sectional

       do k=1,pver
           do i=1,ncol

                condensationSink(:,:) = 0.0_r8  !Sink to the coming "receiver" of any vapour
                !smb++ sectional
                ! initialize condensation sink for sectional
                condensationSink_sec(:,:) = 0.0_r8  !Sink to the coming "receiver" of any vapour
                !smb-- sectional

                !NB: The following is duplicated code, coordinate with koagsub!!
                !Initialize number concentration for this receiver

                !Air density
                rhoAir(i,k) = pmid(i,k)/rair/temperature(i,k)


                numberConcentration(:) = 0.0_r8

                !Go though all modes receiving condensation
                do mode_index_receiver = 0, nmodes

                   !Go through all core species in that mode
                   do tracerIndex = 1, getNumberOfBackgroundTracersInMode(mode_index_receiver)

                      !Find the lifecycle-specie receiving the condensation
                      l_receiver = getTracerIndex(mode_index_receiver, tracerIndex, .true.)

                      !Add up the number concentration of the receiving mode [#/m3]
                      numberConcentration(mode_index_receiver) = numberConcentration(mode_index_receiver) &  !previous value
                                                    + q(i,k,l_receiver)                   &  !kg/kg
                                                    / rhopart(physicsIndex(l_receiver))   &  !m3/kg ==> m3_{aer}/kg_{air}
                                                    * volumeToNumber(mode_index_receiver) &  !#/m3 ==> #/kg_{air}
                                                    * rhoAir(i,k)                                 !kg/m3 ==> #/m3_{air}
                   end do !Lifecycle "core" species in this mode
                enddo
                !smb++ sectional
                ! Calculate number concentration in each bin :
                numberConcentration_sec(i,k, :) = 0.0_r8

                !Go though all bins receiving condensation
                do indBin = 1, secNrBins
                   !Go through all the species in the bin
                   do  indSpec = 1,secNrSpec
                      ind_sec=chemistryIndex(secConstIndex(indSpec, indBin))
                      call sec_numberConc(q(i,k,ind_sec),indSpec, &!SpecNames(tracerIndex), &
                               indBin,rhoAir(i,k), &
                               dummy)        ! dummy holds number concentration
                      ! add number concentration from tracer
                      numberConcentration_sec(i,k, indBin)=numberConcentration_sec(i,k, indBin)+&
                                                                       dummy

                   end do !
                enddo
                !smb-- sectional


                !All modes are condensation receivers
                do cond_vap_idx=1,N_COND_VAP
                   do mode_index_receiver = 0, nmodes

                      !This is the loss rate a gas molecule will see due to aerosol surface area
                      condensationSink(mode_index_receiver,cond_vap_idx)   = &
                              normalizedCondensationSink(mode_index_receiver,cond_vap_idx)  &   ![m3/#/s]
                                      * numberConcentration(mode_index_receiver)                ![#/m3]
                                                                                                !==> [1/s]
                   end do !Loop over receivers
                end do
                !smb++ sectional
                ! condensation sink to sectional bin:
                do indSpec=1,secNrSpec
                   do indBin = 1, secNrBins

                      !This is the loss rate a gas molecule will see due to aerosol surface area
                      condensationSink_sec( indSpec, indBin)   = normalizedCondensationSink_sec( indSpec, indBin)  & ![m3/#/s]  per particle
                                                          * numberConcentration_sec(i,k,indBin)             ![#/m3]
                                                          !==> [1/s]
                   end do !Loop over receivers
                end do
                !smb-- sectional


                !Find concentration after condensation of all
                !condenseable vapours
                ! smb++ sectional edited:
                ! condensation sink fraction that goes to each bin in sectional scheme
                ! set to zero first:
                condensationSinkFraction_sec(i,k,:,:)=0.0_r8
                ! smb-- sectional
                do cond_vap_idx=1,N_COND_VAP

                    !sum of cond. sink for this vapour [1/s]
                    sumCondensationSink(i,k,cond_vap_idx) = sum(condensationSink(:,cond_vap_idx))
                    !smb++ sectional
                    ! Need to add condensation sink to sectional scheme for particles in sectional scheme
                    ! However, not all tracers may contribute:
                    ! assumes same order of gasses in sectional and other (1: H2SO4,2: SOA_LV, 3: SOA_SV)
                    if (cond_vap_idx .LE. secNrSpec) then
                        sumCondensationSink(i,k,cond_vap_idx) = sumCondensationSink(i,k,cond_vap_idx)+&
                                sum(condensationSink_sec( cond_vap_idx,:))
                        ! Keeps track of the fraction of the condensate to the sectional bins for each tracer
                        condensationSinkFraction_sec(i,k,cond_vap_idx,:) = condensationSink_sec(cond_vap_idx, :) &
                                /(sumCondensationSink(i,k,cond_vap_idx)+1.e-30_r8)![frc]
                    end if
                    !smb-- sectional


                !Solve the intermediate (end of timestep) concentration using
                !euler backward solution C_{old} + P *dt - L*C_{new}*dt = C_{new} ==>
                !Cnew -Cold = prod - loss ==>
                intermediateConcentration(i,k,cond_vap_idx) = &
                                     ( q(i,k,cond_vap_map(cond_vap_idx)) + cond_vap_gasprod(i,k,cond_vap_idx)*dt ) &
                                     / (1.0_r8 + sumCondensationSink(i,k,cond_vap_idx)*dt)
                end do

                !Save the fraction of condensation sink for the externally mixed modes
                !(Needed below to find volume shell)
                do cond_vap_idx=1,N_COND_VAP

                   do iDonor = 1,numberOfExternallyMixedModes
                      !Find the mode in question
                      mode_index_donor    = externallyMixedMode(iDonor)

                      !Remember fraction of cond sink for this mode
                      condensationSinkFraction(i,k,iDonor,cond_vap_idx) = &
                            condensationSink(mode_index_donor,cond_vap_idx)    &
                            / sumCondensationSink(i,k,cond_vap_idx)

                      !Remember number concentration in this mode
                      numberConcentrationExtMix(i,k,iDonor) =  &
                             numberConcentration(mode_index_donor)
                   end do
                end do

                !Assume only a fraction of ORG_LV left can contribute to nucleation
                soa_lv_forNucleation(i,k) = lvocfrac*intermediateConcentration(i,k,COND_VAP_ORG_LV) !fraction of soa_lv left that is assumend to have low enough
                                                   !volatility to nucleate.

                modeIndexReceiverCoag = 0
                !Sum coagulation sink for nucleated so4 and soa particles over all receivers of coagulate. Needed for RM's nucleation code
                !OBS - looks like RM's coagulation sink is multiplied by 10^-12??
                do iCoagReceiver = 1, numberOfCoagulationReceiversNPF

                   modeIndexReceiverCoag = receiverModeNPF(iCoagReceiver)

                   coagulationSink(i,k) =   &                                                ![1/s]
                      coagulationSink(i,k) + &                                               ![1/] previous value
                      !++SMB: add npf coagulation
                      !normalizedCoagulationSink(modeIndexReceiverCoag,MODE_IDX_SO4SOA_AIT) & ![m3/#/s]
                      normalizedCoagulationSinkNPF(modeIndexReceiverCoag) & ![m3/#/s]
                                     * numberConcentration(modeIndexReceiverCoag)            !numberConcentration (#/m3)
                      !--SMB: add npf coagulation
                end do    !coagulation sink


           end do !index i
       end do !index k

       !Calculate nucleated masses of so4 and soa (nuclso4, nuclsoa)
       !following RM's parameterization (cka)
       call aeronucl(lchnk,ncol,temperature, pmid, qh20, &
                   intermediateConcentration(:,:,COND_VAP_H2SO4), soa_lv_forNucleation, &
                   coagulationSink, nuclso4, nuclsoa, zm, pblh, &
                   !smb++ sectional
                   nuclrate,nuclrate_pbl_o, formrate, formrate_pbl_o, &
                   orgnucl_o, h2so4nucl_o, grsoa_o, grh2so4_o, dt &
                   ,secMeanD(1) &
                   !smb-- sectional
               )

       !smb++ sectional
       coagnucl_o(:,:)= coagnucl_o(:,:) + coagulationSink(:,:)*dt
       !smb-- sectional
       firstOrderLossRateNucl(:,:,:) = 0.0_r8
       do k=1,pver
          do i=1,ncol

             !First order loss rate (1/s) for nucleation
              !smb++ added check for 0
              if (intermediateConcentration(i,k,COND_VAP_H2SO4) .eq. 0) then
                  firstOrderLossRateNucl(i,k,COND_VAP_H2SO4)=0._r8
              else
                  firstOrderLossRateNucl(i,k,COND_VAP_H2SO4) = nuclSo4(i,k)/intermediateConcentration(i,k,COND_VAP_H2SO4)

              end if
              if (intermediateConcentration(i,k,COND_VAP_ORG_LV) .eq. 0 )then
                  firstOrderLossRateNucl(i,k,COND_VAP_ORG_LV) = 0._r8!nuclSo4(i,k)/(intermediateConcentration(i,k,COND_VAP_ORG_LV))
                else
                  firstOrderLossRateNucl(i,k,COND_VAP_ORG_LV) = nuclSOA(i,k)/intermediateConcentration(i,k,COND_VAP_ORG_LV)
              end if
              !smb-- added check for 0
             !First order loss rate (1/s) for nucleation

             do cond_vap_idx = 1,N_COND_VAP
                !Solve implicitly (again)
                !C_new - C_old =  PROD_{gas} - CS*C_new*dt - LR_{nucl}*C_new =>
                intermediateConcentration(i,k,cond_vap_idx) = &
                               ( q(i,k,cond_vap_map(cond_vap_idx)) + cond_vap_gasprod(i,k,cond_vap_idx)*dt ) &
                               / (1.0_r8 + sumCondensationSink(i,k,cond_vap_idx)*dt + firstOrderLossRateNucl(i,k,cond_vap_idx)*dt)

                !fraction nucleated
                fracNucl(i,k,cond_vap_idx) = firstOrderLossRateNucl(i,k,cond_vap_idx) &
                                     /(firstOrderLossRateNucl(i,k,cond_vap_idx) + sumCondensationSink(i,k,cond_vap_idx))
                !From budget, we get: lost = prod -cnew + cold
                gasLost(i,k,cond_vap_idx) = cond_vap_gasprod(i,k,cond_vap_idx)*dt   & !Produced
                                     + q(i,k,cond_vap_map(cond_vap_idx))            & !cold
                                     - intermediateConcentration(i,k,cond_vap_idx)    !cnew

             end do !cond_vap_idx

             !Add nuceated mass to so4_na mode
             !smb++ sectional don't add to so4_na directly, must go to sectional scheme (later)
             !q(i,k,chemistryIndex(l_so4_na)) =  q(i,k,chemistryIndex(l_so4_na))       &
             !            + gasLost(i,k,COND_VAP_H2SO4)*fracNucl(i,k,COND_VAP_H2SO4)

             !H2SO4 condensate
             do ind_sec=1, secNrBins
                    q(i,k,chemistryIndex(secConstIndex(1,ind_sec))) =                               &
                            q(i,k,chemistryIndex(secConstIndex(1,ind_sec)))                         &
                            + gasLost(i,k,COND_VAP_H2SO4)*(1.0_r8-fracNucl(i,k,COND_VAP_H2SO4))     &
                            *condensationSinkFraction_sec(i,k,COND_VAP_H2SO4, ind_sec)  ! fraction to the particular bin
             end do
             !smb-- sectional
             !H2SO4 condensate
             q(i,k,chemistryIndex(l_so4_a1)) = q(i,k,chemistryIndex(l_so4_a1))         &
                            + gasLost(i,k,COND_VAP_H2SO4)*(1.0_r8-fracNucl(i,k,COND_VAP_H2SO4)) &
                            !smb++ sectional must substract the fraction which goes to the sectional particles:
                            *(1-sum(condensationSinkFraction_sec(i,k,COND_VAP_H2SO4,:)))
                            !smb-- sectional

             !Add nucleated mass to soa_na mode
             !smb++sectional sectional don't add to so4_na directly, must go to sectional scheme (done later)
             !q(i,k,chemistryIndex(l_soa_na)) =  q(i,k,chemistryIndex(l_soa_na))       &
             !            + gasLost(i,k,COND_VAP_ORG_LV)*fracNucl(i,k,COND_VAP_ORG_LV)
             !SOA_LV condensate
             do ind_sec=1, secNrBins
                    q(i,k,chemistryIndex(secConstIndex(2,ind_sec))) =                                   &
                            q(i,k,chemistryIndex(secConstIndex(2,ind_sec)))                             &
                            + gasLost(i,k,COND_VAP_ORG_LV)*(1.0_r8 - fracNucl(i,k,COND_VAP_ORG_LV))     &
                            *condensationSinkFraction_sec(i,k,COND_VAP_ORG_LV, ind_sec) ! fraction to the particular bin
             end do
             !smb--sectional

             !Organic condensate (from both soa_lv and soa_sv) goes to the soaCondensateReceiver tracer (cka)
             q(i,k,chemistryIndex(l_soa_a1)) = q(i,k,chemistryIndex(l_soa_a1))         &
                            + gasLost(i,k,COND_VAP_ORG_SV)                             &           ! "semi volatile" can not nucleate
                            + gasLost(i,k,COND_VAP_ORG_LV)*(1.0_r8-fracNucl(i,k,COND_VAP_ORG_LV)) &  ! part of low volatile which does not nucleate
                            !smb++ sectional must substract the fraction which goes to the sectional particles:
                             *(1-sum(condensationSinkFraction_sec(i,k,COND_VAP_ORG_LV,:)))
                            !smb--sectional

             !condenseable vapours
             q(i,k,chemistryIndex(l_h2so4))  = intermediateConcentration(i,k,COND_VAP_H2SO4)
             q(i,k,chemistryIndex(l_soa_lv)) = intermediateConcentration(i,k,COND_VAP_ORG_LV)
             q(i,k,chemistryIndex(l_soa_sv)) = intermediateConcentration(i,k,COND_VAP_ORG_SV)


             !smb++sectional grow particles in sectional scheme:
             call sec_moveMass(q(i,k,:), numberConcentration_sec(i,k,:), leaveSec(i,k,:), &
                            rhoAir(i,k), 2._r8*originalNumberMedianRadius(MODE_IDX_SO4SOA_AIT), split_dt)
             ! Add nucleated mass to first bin of sectional scheme:
             q(i,k,chemistryIndex(secConstIndex(1,1))) =  q(i,k,chemistryIndex(secConstIndex(1,1)))       &
                         + gasLost(i,k,COND_VAP_H2SO4)*fracNucl(i,k,COND_VAP_H2SO4)
             q(i,k,chemistryIndex(secConstIndex(2,1))) =  q(i,k,chemistryIndex(secConstIndex(2,1)))       &
                         + gasLost(i,k,COND_VAP_ORG_LV)*fracNucl(i,k,COND_VAP_ORG_LV)

            ! Add mass from sectional scheme to so4_na and soa_na:
             q(i,k,chemistryIndex(l_so4_na)) = q(i,k,chemistryIndex(l_so4_na))         &
                                    +leaveSec(i,k, 1)
             q(i,k,chemistryIndex(l_soa_na)) = q(i,k,chemistryIndex(l_soa_na))         &
                                    +leaveSec(i,k,2)
             !smb--sectional


             !Condensation transfers mass from externally mixed to internally mixed modes
             do iDonor = 1,numberOfExternallyMixedModes

                !Find the mode in question
                mode_index_donor    = externallyMixedMode(iDonor)

                if(getNumberOfTracersInMode(mode_index_donor) .eq. 0)then
                   cycle
                end if

                volume_shell = 0.0_r8
                do cond_vap_idx = 1, N_COND_VAP

                   !Add up volume shell for this
                   !condenseable vapour
                   volume_shell = volume_shell                                               &
                         + condensationSinkFraction(i,k,iDonor,cond_vap_idx)                 & ![frc]
                         * gasLost(i,k,cond_vap_idx)*(1.0_r8-fracNucl(i,k,cond_vap_idx))     & ![kg/kg]
                         * invRhoPart(physicsIndex(cond_vap_map(cond_vap_idx)))              & !*[m3/kg] ==> [m3/kg_{air}
                         * rhoAir(i,k)                                                         !*[kg/m3] ==> m3/m3

                end do

                area_core=numberConcentrationExtMix(i,k,iDonor)*numberToSurface(mode_index_donor)   !#/m3 * m2/# ==> m2/m3
                vol_monolayer=area_core*dr_so4_monolayers_age

                ! Small fraction retained to avoid numerical irregularities
                frac_transfer=min((volume_shell/vol_monolayer),0.999_r8)

                !How many tracers exist in donor mode?
                !The "donor" is the externally mixed mode which will soon
                !become internally mixed. The externally mixed is donating mass
                !and the internally mixed is receiving...
                do tracerIndex = 1, getNumberOfTracersInMode(mode_index_donor)

                   !Indexes here are in "chemistry space"
                   l_donor    = getTracerIndex(mode_index_donor, tracerIndex,.true.)
                   l_receiver = lifeCycleReceiver(l_donor)

                   if( l_receiver .le. 0)then
                      stop !something wrong
                   endif

                   !Transfer from donor to receiver takes into account
                   !fraction transferred
                   totalLoss(i,k,l_donor) = frac_transfer*q(i,k,l_donor)
                   q(i,k,l_donor) = q(i,k,l_donor) - totalLoss(i,k,l_donor)
                   q(i,k,l_receiver) = q(i,k,l_receiver) + totalLoss(i,k,l_donor)
                end do !tracers in mode
             end do    !loop over receivers
          end do !physical index k
       end do    !physical index i

       !Output for diagnostics
       call phys_getopts(history_aerosol_out = history_aerosol)

       if(history_aerosol)then
          coltend(:ncol,:) = 0.0_r8
          do i=1,gas_pcnst
             !Check if species contributes to condensation
             if(lifeCycleReceiver(i) .gt. 0)then
               !Loss from the donor specie
               tracer_coltend(:ncol) = sum(totalLoss(:ncol, :,i)*pdel(:ncol,:),2)/gravit/dt
               coltend(:ncol,i) = coltend(:ncol,i) - tracer_coltend(:ncol) !negative (loss for donor)
               coltend(:ncol,lifeCycleReceiver(i)) = coltend(:ncol,lifeCycleReceiver(i)) + tracer_coltend(:ncol)
             endif
          end do


          !smb++ sectional
          ! Remove so4_n ---> directly into so4_na
          coltend(:ncol,chemistryIndex(secConstIndex(1,1))) = coltend(:ncol,chemistryIndex(secConstIndex(1,1))) + &
          !smb-- sectional
                                                 sum(                                         &
                                                    gasLost(:ncol,:,COND_VAP_H2SO4)           &
                                                    *fracNucl(:ncol,:,COND_VAP_H2SO4)*pdel(:ncol,:) , 2 &
                                                    )/gravit/dt

          ! Remove so4_n ---> directly into so4_na
          !smb++ put in how much leaves sectional scheme:
          coltend(:ncol,chemistryIndex(l_so4_na)) = coltend(:ncol,chemistryIndex(l_so4_na)) + &
                                                 sum(                                         &
                                                    leaveSec(:ncol,:,COND_VAP_H2SO4)           &
                                                    *pdel(:ncol,:) , 2 &
                                                    )/gravit/dt
          !smb--sectional
          !Take into account H2SO4 (gas) condensed in budget
          coltend(:ncol,chemistryIndex(l_so4_a1)) = coltend(:ncol,chemistryIndex(l_so4_a1)) + &
                                                 sum(                                         &
                                                    gasLost(:ncol,:,COND_VAP_H2SO4)           &
                                                    !smb++sectional subtract fraction to sectional scheme
                                                    *(1-sum(condensationSinkFraction_sec(:ncol,:,COND_VAP_H2SO4,:),3)) &
                                                    !smb--sectional
                                                    *(1.0_r8 - fracNucl(:ncol,:,COND_VAP_H2SO4))*pdel(:ncol,:) , 2 &
                                                    )/gravit/dt

          !Take into account soa_lv (gas) nucleated in budget
          !smb++ sectional
          coltend(:ncol,chemistryIndex(secConstIndex(2,1))) = coltend(:ncol,chemistryIndex(secConstIndex(2,1))) + &
          !smb-- sectional
                                                 sum(                                         &
                                                    gasLost(:ncol,:,COND_VAP_ORG_LV)              &
                                                    *fracNucl(:ncol,:,COND_VAP_ORG_LV)*pdel(:ncol,:) , 2 &
                                                    )/gravit/dt

          !Take into account soa_lv (gas) nucleated in budget
          !smb++ sectional  put in how much leaves sectional scheme:
          coltend(:ncol,chemistryIndex(l_soa_na)) = coltend(:ncol,chemistryIndex(l_soa_na)) + &
                                                 sum(                                         &
                                                    leaveSec(:ncol,:,COND_VAP_ORG_LV)              &
                                                    *pdel(:ncol,:) , 2 &
                                                    )/gravit/dt
          !smb-- sectional
          !Take into account soa gas condensed in the budget (both LV and SV)
          coltend(:ncol,chemistryIndex(l_soa_a1)) = coltend(:ncol,chemistryIndex(l_soa_a1)) + &
                                                 sum(                                         &
                                                    gasLost(:ncol,:,COND_VAP_ORG_LV)           &
                                                    !smb++sectional subtract fraction to sectional scheme
                                                    *(1-sum(condensationSinkFraction_sec(:ncol,:,COND_VAP_ORG_LV,:),3)) &
                                                    !smb--sectional
                                                    *(1.0_r8 - fracNucl(:ncol,:,COND_VAP_ORG_LV))*pdel(:ncol,:) , 2 &
                                                    )/gravit/dt                        &
                                                    +                                  &
                                                 sum(                                         &
                                                    gasLost(:ncol,:,COND_VAP_ORG_SV)*pdel(:ncol,:) , 2 &
                                                    )/gravit/dt
          !smb++ sectional: putt in how condenses on sectional:
          do indBin=1,secNrBins
               coltend(:ncol,chemistryIndex(secConstIndex(1,indBin))) = coltend(:ncol, chemistryIndex(secConstIndex(1,indBin)))+ &
                                                    sum(                    &
                                                    gasLost(:ncol, :, COND_VAP_H2SO4)    &
                                                    *(condensationSinkFraction_sec(:ncol,:,COND_VAP_H2SO4,indBin)) &
                                                    *(1.0_r8 - fracNucl(:ncol, :, COND_VAP_H2SO4))*pdel(:ncol,:),2) &
                                                    /gravit/dt
          end do
          do indBin=1,secNrBins
               coltend(:ncol,chemistryIndex(secConstIndex(2,indBin))) = coltend(:ncol, chemistryIndex(secConstIndex(2,indBin)))+ &
                                                    sum(                    &
                                                    gasLost(:ncol, :, COND_VAP_ORG_LV)    &
                                                    *(condensationSinkFraction_sec(:ncol,:,COND_VAP_ORG_LV,indBin)) &
                                                    *(1.0_r8 - fracNucl(:ncol, :, COND_VAP_ORG_LV))*pdel(:ncol,:),2) &
                                                    /gravit/dt
          end do

          coltend_o(:,:)=coltend(:,:)
          !smb-- sectional
          !smb++ sectional remove write to file, done in super method
          ! do i=1,gas_pcnst
          !    if(lifeCycleReceiver(i) .gt. 0 )then
          !       long_name= trim(solsym(i))//"condTend"
          !       call outfld(long_name, coltend(:ncol,i), pcols, lchnk)
          !       long_name= trim(solsym(lifeCycleReceiver(i)))//"condTend"
          !       call outfld(long_name, coltend(:ncol,lifeCycleReceiver(i)),pcols,lchnk)
          !    end if
          ! end do
          ! long_name=trim(solsym(chemistryIndex(l_so4_a1)))//"condTend"
          ! call outfld(long_name, coltend(:ncol,chemistryIndex(l_so4_a1)),pcols,lchnk)
          ! long_name=trim(solsym(chemistryIndex(l_soa_a1)))//"condTend"
          ! call outfld(long_name, coltend(:ncol,chemistryIndex(l_soa_a1)),pcols,lchnk)
          ! long_name=trim(solsym(chemistryIndex(l_so4_na)))//"condTend"
          ! call outfld(long_name, coltend(:ncol,chemistryIndex(l_so4_na)),pcols,lchnk)
          ! long_name=trim(solsym(chemistryIndex(l_soa_na)))//"condTend"
          ! call outfld(long_name, coltend(:ncol,chemistryIndex(l_soa_na)),pcols,lchnk)
          !smb-- sectional

       endif


       return
   end subroutine condtend_sub


end module condtend