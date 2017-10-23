module Reaction_Sandbox_S2o4_fe3_class

! # Id: reaction_sandbox_s2o4_fe33, Tue 08 Mar 2016 04:43:33 PM MST pandeys #
!       - update s2o4, so3, hcro4, hco3, h, cr(oh)3(s), fe(oh)3(s), siderite
!       - full reactant inhibition, no product inhibition

! 1. Change all references to "S2o4_fe3" as desired to rename the module and
!    and subroutines within the module.

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

! 2. Add module variables here.  Note that one must use the PETSc data types
!    PetscInt, PetscReal, PetscBool to declare variables of type integer
!    float/real*8, and logical respectively.  E.g.
!
! PetscReal, parameter :: formula_weight_of_water = 18.01534d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_s2o4_fe3_type
! 3. Add variables/arrays associated with new reaction.
    character(len=MAXWORDLENGTH) :: name_spec_s2o4 ! S2O4--
    character(len=MAXWORDLENGTH) :: name_spec_h ! H+
    character(len=MAXWORDLENGTH) :: name_spec_so3 ! SO3--
    character(len=MAXWORDLENGTH) :: name_mnrl_fe3 ! Fe(OH)3(s)
    character(len=MAXWORDLENGTH) :: name_bound_fe2_slow ! bound_Fe++ SLOW
    character(len=MAXWORDLENGTH) :: name_bound_fe2_fast ! bound_Fe++ FAST
    PetscInt :: id_spec_s2o4, id_spec_h, id_spec_so3
    PetscInt :: id_mnrl_fe3, id_bound_fe2_slow, id_bound_fe2_fast
    PetscReal :: rate_constant, ssa, rock_density, sitefrac, eps
  contains
    procedure, public :: ReadInput => S2o4_fe3Read
    procedure, public :: Setup => S2o4_fe3Setup
    procedure, public :: Evaluate => S2o4_fe3React
    procedure, public :: UpdateKineticState => S2o4_fe3KineticState
    procedure, public :: Destroy => S2o4_fe3Destroy
  end type reaction_sandbox_s2o4_fe3_type

  public :: S2o4_fe3Create

contains

! ************************************************************************** !

function S2o4_fe3Create()
  !
  ! Allocates s2o4_fe3 reaction object.
  !
  ! Author: Sachin Pandey 
  ! Date: 09/08/16 
  !

  implicit none

  class(reaction_sandbox_s2o4_fe3_type), pointer :: S2o4_fe3Create

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(S2o4_fe3Create)
  S2o4_fe3Create%name_spec_s2o4 = ''
  S2o4_fe3Create%name_spec_h = ''
  S2o4_fe3Create%name_spec_so3 = ''
  S2o4_fe3Create%name_mnrl_fe3 = ''
  S2o4_fe3Create%name_bound_fe2_slow = ''
  S2o4_fe3Create%name_bound_fe2_fast = ''
  S2o4_fe3Create%id_spec_s2o4 = 0
  S2o4_fe3Create%id_spec_h = 0
  S2o4_fe3Create%id_spec_so3 = 0
  S2o4_fe3Create%id_mnrl_fe3 = 0
  S2o4_fe3Create%id_bound_fe2_slow = 0
  S2o4_fe3Create%id_bound_fe2_fast = 0
  S2o4_fe3Create%rate_constant = 0.d0
  S2o4_fe3Create%ssa = 0.d0
  S2o4_fe3Create%rock_density = 0.d0
  S2o4_fe3Create%sitefrac = 0.d0
  S2o4_fe3Create%eps = 0.d0

  nullify(S2o4_fe3Create%next)

end function S2o4_fe3Create

! ************************************************************************** !

subroutine S2o4_fe3Read(this,input,option)
  !
  ! Reads input deck for s2o4_fe3 reaction parameters (if any)
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_s2o4_fe3_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,S2O4_FE3')
    call StringToUpper(word)

    select case(trim(word))

      ! S2o4_fe3 Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     S2O4_FE3
      !       S2O4_FE3_INTEGER 1
      !       S2O4_FE3_INTEGER_ARRAY 2 3 4
      !     END
      !   : end user defined input
      !   END
      !   ...
      ! END

! 5. Add case statement for reading variables.  E.g.
! 6. Read the variable
! 7. Inform the user of any errors if not read correctly.
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'RATE_CONSTANT', &
                           'CHEMISTRY,REACTION_SANDBOX,S2O4_FE3')
      case('SSA')
        call InputReadDouble(input,option,this%ssa)
        call InputErrorMsg(input,option,'SSA', &
                           'CHEMISTRY,REACTION_SANDBOX,S2O4_FE3')
      case('ROCK_DENSITY')
        call InputReadDouble(input,option,this%rock_density)
        call InputErrorMsg(input,option,'ROCK_DENSITY', &
                           'CHEMISTRY,REACTION_SANDBOX,S2O4_FE3')
      case('FRACTION')
        call InputReadDouble(input,option,this%sitefrac)
        call InputErrorMsg(input,option,'FRACTION', &
                           'CHEMISTRY,REACTION_SANDBOX,S2O4_FE3')
      case('EPS')
        call InputReadDouble(input,option,this%eps)
        call InputErrorMsg(input,option,'EPS', &
                           'CHEMISTRY,REACTION_SANDBOX,S2O4_FE3')
! 8. Inform the user if the keyword is not recognized
      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,S2O4_FE3',option)
    end select
  enddo

end subroutine S2o4_fe3Read

! ************************************************************************** !

subroutine S2o4_fe3Setup(this,reaction,option)
  !
  ! Sets up the s2o4_fe3 reaction either with parameters either
  ! read from the input deck or hardwired.
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_s2o4_fe3_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize
  this%name_spec_s2o4 = 'S2O4--'
  this%name_spec_h    = 'H+'
  this%name_spec_so3  = 'SO3--'
  this%name_mnrl_fe3  = 'Fe(OH)3(s)'
  this%name_bound_fe2_slow  = 'slow_Fe++'
  this%name_bound_fe2_fast  = 'fast_Fe++'
  
  ! Aqueous species
  this%id_spec_s2o4 = &
    GetPrimarySpeciesIDFromName(this%name_spec_s2o4,reaction,option)
  this%id_spec_h = &
      GetPrimarySpeciesIDFromName(this%name_spec_h,reaction,option)
  this%id_spec_so3 = &
      GetPrimarySpeciesIDFromName(this%name_spec_so3,reaction,option)

  ! Mineral species
  this%id_mnrl_fe3 = &
    GetMineralIDFromName(this%name_mnrl_fe3,reaction%mineral,option)

  ! Immobile species
  this%id_bound_fe2_slow = &
    GetImmobileSpeciesIDFromName(this%name_bound_fe2_slow,reaction%immobile,option)
  this%id_bound_fe2_fast = &
    GetImmobileSpeciesIDFromName(this%name_bound_fe2_fast,reaction%immobile,option)


end subroutine S2o4_fe3Setup

! ************************************************************************** !

subroutine S2o4_fe3React(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_s2o4_fe3_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water, FeII, FeIII, rate
  PetscReal :: vf_feoh3, mv_feoh3, mw_feoh3
  PetscInt :: id_bound_fe2_slow_offset, id_bound_fe2_fast_offset

  ! Info for Fe(oh)3(s) and bound_Fe++
  vf_feoh3 = rt_auxvar%mnrl_volfrac(this%id_mnrl_fe3) ! m^3/m^3_bulk
  mv_feoh3 = reaction%mineral%kinmnrl_molar_vol(this%id_mnrl_fe3) ! m^3/mol
  mw_feoh3 = reaction%mineral%kinmnrl_molar_wt(this%id_mnrl_fe3) ! m^3/mol
  id_bound_fe2_slow_offset = reaction%offset_immobile + this%id_bound_fe2_slow
  id_bound_fe2_fast_offset = reaction%offset_immobile + this%id_bound_fe2_fast

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! L_water from m^3_water

  ! UNITS: g_fe(oh)3/g_sed from from m^3_mnrl/m^3_bulk
  FeIII = vf_feoh3/mv_feoh3/this%rock_density/1.d3*mw_feoh3

  ! mole/l-s
  rate = this%rate_constant*this%ssa* &
             rt_auxvar%pri_molal(this%id_spec_s2o4)* &
             FeIII
 
  if (rate < this%eps) then
    rate = 0.d0
  endif

  ! for debug
#if 0
   print *, 's2o4', rt_auxvar%total(this%id_spec_s2o4,iphase)
   print *, 'FeIII', FeIII
   print *, 'vf_feoh3', vf_feoh3
   print *, 'mv_feoh3', mv_feoh3
   print *, 'mw_feoh3', mw_feoh3
   print *, 'rock_density', this%rock_density
#endif

  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second  
  Residual(this%id_spec_s2o4) = Residual(this%id_spec_s2o4) - (-1.0*rate) * L_water
  Residual(this%id_spec_h) = Residual(this%id_spec_h) - (-2.0*rate) * L_water
  Residual(this%id_spec_so3) = Residual(this%id_spec_so3) - (2.0*rate) * L_water
  Residual(id_bound_fe2_slow_offset) = Residual(id_bound_fe2_slow_offset) - (2.0*rate*L_water) * (1-this%sitefrac)
  Residual(id_bound_fe2_fast_offset) = Residual(id_bound_fe2_fast_offset) - (2.0*rate*L_water) * this%sitefrac

  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

     option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
                        'due to assumptions in S2o4_fe3'
     call printErrMsg(option)

  endif

end subroutine S2o4_fe3React

! ************************************************************************** !

subroutine S2o4_fe3KineticState(this,rt_auxvar,global_auxvar, &
                                  material_auxvar,reaction,option)
  !
  ! For update mineral volume fractions
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_s2o4_fe3_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  ! the following arrays must be declared after reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: FeII, FeIII, rate, cnv_L2bulk
  PetscReal :: vf_feoh3, mv_feoh3, mw_feoh3
  PetscReal :: delta_volfrac

  ! Unit of the residual must be in moles/second
  ! Unit of mnr_rate should be mol/sec/m^3 bulk
  ! 1.d3 converts m^3 water -> L water
  vf_feoh3 = rt_auxvar%mnrl_volfrac(this%id_mnrl_fe3) ! m^3/m^3_bulk
  mv_feoh3 = reaction%mineral%kinmnrl_molar_vol(this%id_mnrl_fe3) ! m^3/mol
  mw_feoh3 = reaction%mineral%kinmnrl_molar_wt(this%id_mnrl_fe3) ! m^3/mol

  ! L_water/m^3_bulk
  cnv_L2bulk = material_auxvar%porosity*global_auxvar%sat(iphase)*1.d3

  ! UNITS: m^2_reactant/g_sed from from m^3_mnrl/m^3_bulk
  FeIII = vf_feoh3/mv_feoh3/this%rock_density/1.d3*mw_feoh3

  ! mole/l-s
  rate = this%rate_constant*this%ssa* &
             rt_auxvar%pri_molal(this%id_spec_s2o4)* &
             FeIII

  if (rate < this%eps) then
    rate = 0.d0
  endif

  ! rate = mol/m^3/sec
  ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) * mol_vol (m^3 mnrl/mol mnrl)
  ! Update Fe(OH)3
  delta_volfrac = (-2.0*rate)*cnv_L2bulk* &
                  reaction%mineral%kinmnrl_molar_vol(this%id_mnrl_fe3)* &
                  option%tran_dt
  rt_auxvar%mnrl_volfrac(this%id_mnrl_fe3) = &
    rt_auxvar%mnrl_volfrac(this%id_mnrl_fe3) + delta_volfrac

end subroutine S2o4_fe3KineticState

! ************************************************************************** !

subroutine S2o4_fe3Destroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  implicit none

  class(reaction_sandbox_s2o4_fe3_type) :: this

! 12. Add code to deallocate contents of the s2o4_fe3 object

end subroutine S2o4_fe3Destroy

end module Reaction_Sandbox_S2o4_fe3_class
