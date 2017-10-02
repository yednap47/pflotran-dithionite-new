module Reaction_Sandbox_Fe2_o2_class

! # Id: reaction_sandbox_fe2_o2, Fri 05 May 2016 04:43:33 PM MST pandeys #

! 1. Change all references to "Fe2_o2" as desired to rename the module and
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
    extends(reaction_sandbox_base_type) :: reaction_sandbox_fe2_o2_type
! 3. Add variables/arrays associated with new reaction.
    character(len=MAXWORDLENGTH) :: name_spec_h ! H+
    character(len=MAXWORDLENGTH) :: name_spec_o2 ! O2(aq)
    character(len=MAXWORDLENGTH) :: name_spec_fe3 ! Fe+++
    character(len=MAXWORDLENGTH) :: name_bound_fe2_slow ! bound_Fe++
    character(len=MAXWORDLENGTH) :: name_bound_fe2_fast ! bound_Fe++
    PetscInt :: id_spec_h, id_spec_o2, id_spec_fe3
    PetscInt :: id_bound_fe2_slow, id_bound_fe2_fast
    PetscReal :: rate_constant_slow, rate_constant_fast, rock_density, eps
  contains
    procedure, public :: ReadInput => Fe2_o2Read
    procedure, public :: Setup => Fe2_o2Setup
    procedure, public :: Evaluate => Fe2_o2React
    procedure, public :: UpdateKineticState => Fe2_o2KineticState
    procedure, public :: Destroy => Fe2_o2Destroy
  end type reaction_sandbox_fe2_o2_type

  public :: Fe2_o2Create

contains

! ************************************************************************** !

function Fe2_o2Create()
  !
  ! Allocates fe2_o2 reaction object.
  !
  ! Author: Sachin Pandey 
  ! Date: 05/05/17 
  !

  implicit none

  class(reaction_sandbox_fe2_o2_type), pointer :: Fe2_o2Create

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(Fe2_o2Create)
  Fe2_o2Create%name_spec_h = ''
  Fe2_o2Create%name_spec_o2 = ''
  Fe2_o2Create%name_spec_fe3 = ''
  Fe2_o2Create%name_bound_fe2_slow = ''
  Fe2_o2Create%name_bound_fe2_fast = ''
  Fe2_o2Create%id_spec_h = 0
  Fe2_o2Create%id_spec_o2 = 0
  Fe2_o2Create%id_spec_fe3 = 0
  Fe2_o2Create%id_bound_fe2_slow = 0
  Fe2_o2Create%id_bound_fe2_fast = 0
  Fe2_o2Create%rate_constant_slow = 0.d0
  Fe2_o2Create%rate_constant_fast = 0.d0
  Fe2_o2Create%rock_density = 0.d0
  Fe2_o2Create%eps = 0.d0
  nullify(Fe2_o2Create%next)

end function Fe2_o2Create

! ************************************************************************** !

subroutine Fe2_o2Read(this,input,option)
  !
  ! Reads input deck for fe2_o2 reaction parameters (if any)
  !
  ! Author: Sachin Pandey
  ! Date: 05/05/17
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_fe2_o2_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,FE2_O2')
    call StringToUpper(word)

    select case(trim(word))

      ! Fe2_o2 Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     FE2_O2
      !       FE2_O2_INTEGER 1
      !       FE2_O2_INTEGER_ARRAY 2 3 4
      !     END
      !   : end user defined input
      !   END
      !   ...
      ! END

! 5. Add case statement for reading variables.  E.g.
! 6. Read the variable
! 7. Inform the user of any errors if not read correctly.
     case('RATE_CONSTANT_SLOW')
       call InputReadDouble(input,option,this%rate_constant_slow)
       call InputErrorMsg(input,option,'RATE_CONSTANT_SLOW', &
                           'CHEMISTRY,REACTION_SANDBOX,FE2_O2')
     case('RATE_CONSTANT_FAST')
       call InputReadDouble(input,option,this%rate_constant_fast)
       call InputErrorMsg(input,option,'RATE_CONSTANT_FAST', &
                       'CHEMISTRY,REACTION_SANDBOX,FE2_O2')
     case('ROCK_DENSITY')
       call InputReadDouble(input,option,this%rock_density)
       call InputErrorMsg(input,option,'ROCK_DENSITY', &
                           'CHEMISTRY,REACTION_SANDBOX,FE2_O2')
     case('EPS')
       call InputReadDouble(input,option,this%eps)
       call InputErrorMsg(input,option,'EPS', &
                           'CHEMISTRY,REACTION_SANDBOX,FE2_O2')
! 8. Inform the user if the keyword is not recognized
      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,FE2_O2',option)
    end select
  enddo

end subroutine Fe2_o2Read

! ************************************************************************** !

subroutine Fe2_o2Setup(this,reaction,option)
  !
  ! Sets up the fe2_o2 reaction either with parameters either
  ! read from the input deck or hardwired.
  !
  ! Author: Sachin Pandey
  ! Date: 05/05/17
  !

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_fe2_o2_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize
  this%name_spec_h = 'H+'
  this%name_spec_o2 = 'O2(aq)'
  this%name_spec_fe3 = 'Fe+++'
  this%name_bound_fe2_slow = 'slow_Fe++'
  this%name_bound_fe2_fast = 'fast_Fe++'

  this%id_spec_h = &
      GetPrimarySpeciesIDFromName(this%name_spec_h,reaction,option)
  this%id_spec_o2 = &
      GetPrimarySpeciesIDFromName(this%name_spec_o2,reaction,option)
  this%id_spec_fe3 = &
      GetPrimarySpeciesIDFromName(this%name_spec_fe3,reaction,option)

  ! Immobile species
  this%id_bound_fe2_slow = &
    GetImmobileSpeciesIDFromName(this%name_bound_fe2_slow,reaction%immobile,option)
  this%id_bound_fe2_fast = &
    GetImmobileSpeciesIDFromName(this%name_bound_fe2_fast,reaction%immobile,option)

end subroutine Fe2_o2Setup

! ************************************************************************** !

subroutine Fe2_o2React(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: Sachin Pandey
  ! Date: 05/05/17
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_fe2_o2_type) :: this
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
  PetscReal :: L_water, FeII_slow, FeII_fast, rate_slow, rate_fast
  PetscInt :: id_bound_fe2_slow_offset, id_bound_fe2_fast_offset

  ! Info for bound_Fe++
  id_bound_fe2_slow_offset = reaction%offset_immobile + this%id_bound_fe2_slow
  id_bound_fe2_fast_offset = reaction%offset_immobile + this%id_bound_fe2_fast

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! L_water from m^3_water

  ! UNITS: mol_reactant/g_sed from from mol_reactant/m^3_bulk
  FeII_slow = rt_auxvar%immobile(this%id_bound_fe2_slow)/(this%rock_density*1000)
  FeII_fast = rt_auxvar%immobile(this%id_bound_fe2_fast)/(this%rock_density*1000)

  ! mole/l-s
  rate_slow = this%rate_constant_slow*(FeII_slow)*rt_auxvar%pri_molal(this%id_spec_o2)
  rate_fast = this%rate_constant_fast*(FeII_fast)*rt_auxvar%pri_molal(this%id_spec_o2)

  if (rate_fast < this%eps) then
    rate_fast = 0.d0
  endif

  if (rate_slow < this%eps) then
    rate_slow = 0.d0
  endif

  ! for debug
#if 0
   print *, 'FeII', FeII
   print *, 'O2(aq)', rt_auxvar%pri_molal(this%id_spec_o2)
   print *, 'rate', rate
#endif
 
  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second  
  Residual(this%id_spec_o2) = Residual(this%id_spec_o2) - (-0.25*rate_slow-0.25*rate_fast) * L_water
  Residual(this%id_spec_h) = Residual(this%id_spec_h) - (-1.0*rate_slow-1.0*rate_fast) * L_water
  Residual(this%id_spec_fe3) = Residual(this%id_spec_fe3) - (1.0*rate_slow+1.0*rate_fast) * L_water
  Residual(id_bound_fe2_slow_offset) = Residual(id_bound_fe2_slow_offset) - (-1.0*rate_slow) * L_water
  Residual(id_bound_fe2_fast_offset) = Residual(id_bound_fe2_fast_offset) - (-1.0*rate_fast) * L_water

  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

     option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
                        'due to assumptions in Fe2_o2'
     call printErrMsg(option)

  endif

end subroutine Fe2_o2React

! ************************************************************************** !

subroutine Fe2_o2KineticState(this,rt_auxvar,global_auxvar, &
                                  material_auxvar,reaction,option)
  !
  ! For update mineral volume fractions
  !
  ! Author: Sachin Pandey
  ! Date: 05/05/17
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_fe2_o2_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  ! the following arrays must be declared after reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1

end subroutine Fe2_o2KineticState

! ************************************************************************** !

subroutine Fe2_o2Destroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Sachin Pandey
  ! Date: 05/05/17
  !

  implicit none

  class(reaction_sandbox_fe2_o2_type) :: this

! 12. Add code to deallocate contents of the fe2_o2 object

end subroutine Fe2_o2Destroy

end module Reaction_Sandbox_Fe2_o2_class
