module Reaction_Sandbox_S2o4_disp_class

! Id: Reaction_Sandbox_S2o4_disp, Mon 24 Apr 2017 12:11:58 PM MDT pandeys !
! Created by Sachin Pandey, LANL
! Description: Use first order decay for S2O4 in contact with sediments
!   From Amonette (1994):
!     -dC/dt = k1C + k2CsurfC, k2 >> k1
!     No information to constraint Csurf, simplify to Istok (1999):
!     -dC/dt = k2*C (i.e. first order reaction)
!
!   Fruchter (2000):
!     2S2O4-- + H2O -> 2SO3-- + S2O3-- + H+     
!------------------------------------------------------------------------------

! 1. Change all references to "S2o4_disp" as desired to rename the module and
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
    extends(reaction_sandbox_base_type) :: reaction_sandbox_s2o4_disp_type
! 3. Add variables/arrays associated with new reaction.
    character(len=MAXWORDLENGTH) :: name_spec_s2o4 ! S2O4--
    character(len=MAXWORDLENGTH) :: name_spec_so3 ! SO3--
    character(len=MAXWORDLENGTH) :: name_spec_s2o3  ! S2O3--
    character(len=MAXWORDLENGTH) :: name_spec_h ! H+

    PetscInt :: id_spec_s2o4, id_spec_so3, id_spec_s2o3, id_spec_h
    PetscReal :: rate_constant, eps
  contains
    procedure, public :: ReadInput => S2o4_dispRead
    procedure, public :: Setup => S2o4_dispSetup
    procedure, public :: Evaluate => S2o4_dispReact
    procedure, public :: UpdateKineticState => S2o4_dispKineticState
    procedure, public :: Destroy => S2o4_dispDestroy
  end type reaction_sandbox_s2o4_disp_type

  public :: S2o4_dispCreate

contains

! ************************************************************************** !

function S2o4_dispCreate()
  !
  ! Allocates s2o4_disp reaction object.
  !
  ! Author: Sachin Pandey 
  ! Date: 04/24/17 
  !

  implicit none

  class(reaction_sandbox_s2o4_disp_type), pointer :: S2o4_dispCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(S2o4_dispCreate)
  S2o4_dispCreate%name_spec_s2o4 = ''
  S2o4_dispCreate%name_spec_so3 = ''
  S2o4_dispCreate%name_spec_s2o3 = ''
  S2o4_dispCreate%name_spec_h = ''
  S2o4_dispCreate%id_spec_s2o4 = 0
  S2o4_dispCreate%id_spec_so3 = 0
  S2o4_dispCreate%id_spec_s2o3 = 0
  S2o4_dispCreate%id_spec_h = 0
  S2o4_dispCreate%rate_constant = 0.d0
  S2o4_dispCreate%eps = 0.d0

  nullify(S2o4_dispCreate%next)

end function S2o4_dispCreate

! ************************************************************************** !

subroutine S2o4_dispRead(this,input,option)
  !
  ! Reads input deck for s2o4_disp reaction parameters (if any)
  !
  ! Author: Sachin Pandey
  ! Date: 04/24/17
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_s2o4_disp_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,S2O4_DISP')
    call StringToUpper(word)

    select case(trim(word))

      ! S2o4_disp Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     S2O4_DISP
      !       S2O4_DISP_INTEGER 1
      !       S2O4_DISP_INTEGER_ARRAY 2 3 4
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
                           'CHEMISTRY,REACTION_SANDBOX,S2O4_DISP')
      case('EPS')
        call InputReadDouble(input,option,this%eps)
        call InputErrorMsg(input,option,'EPS', &
                           'CHEMISTRY,REACTION_SANDBOX,S2O4_DISP')
! 8. Inform the user if the keyword is not recognized
      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,S2O4_DISP',option)
    end select
  enddo

end subroutine S2o4_dispRead

! ************************************************************************** !

subroutine S2o4_dispSetup(this,reaction,option)
  !
  ! Sets up the s2o4_disp reaction either with parameters either
  ! read from the input deck or hardwired.
  !
  ! Author: Sachin Pandey
  ! Date: 04/24/17
  !

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_s2o4_disp_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize
  this%name_spec_s2o4 = 'S2O4--'
  this%name_spec_so3 = 'SO3--'
  this%name_spec_s2o3  = 'S2O3--'
  this%name_spec_h = 'H+'

  this%id_spec_s2o4 = &
    GetPrimarySpeciesIDFromName(this%name_spec_s2o4,reaction,option)
  this%id_spec_so3 = &
    GetPrimarySpeciesIDFromName(this%name_spec_so3,reaction,option)
  this%id_spec_s2o3 = &
    GetPrimarySpeciesIDFromName(this%name_spec_s2o3,reaction,option)
  this%id_spec_h = &
    GetPrimarySpeciesIDFromName(this%name_spec_h,reaction,option)

end subroutine S2o4_dispSetup

! ************************************************************************** !

subroutine S2o4_dispReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: Sachin Pandey
  ! Date: 04/24/17
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_s2o4_disp_type) :: this
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
  PetscReal :: L_water, rate

  ! Unit of the residual must be in moles/second

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! L_water from m^3_water

  rate  = this%rate_constant*rt_auxvar%pri_molal(this%id_spec_s2o4)! mole/l-s

  if (rate < this%eps) then
    rate = 0.d0
  endif

  ! Residuals
  Residual(this%id_spec_s2o4) = Residual(this%id_spec_s2o4) - (-1.0*rate) * L_water
  Residual(this%id_spec_so3) = Residual(this%id_spec_so3) - (1.0*rate) * L_water
  Residual(this%id_spec_s2o3) = Residual(this%id_spec_s2o3) - (0.5*rate) * L_water
  Residual(this%id_spec_h) = Residual(this%id_spec_h) - (1.0*rate) * L_water

  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

     option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
                        'due to assumptions in S2o4_disp'
     call printErrMsg(option)

  endif

end subroutine S2o4_dispReact

! ************************************************************************** !

subroutine S2o4_dispKineticState(this,rt_auxvar,global_auxvar, &
                                  material_auxvar,reaction,option)
  !
  ! For update mineral volume fractions
  !
  ! Author: Sachin Pandey
  ! Date: 04/24/17
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_s2o4_disp_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  ! the following arrays must be declared after reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1

end subroutine S2o4_dispKineticState

! ************************************************************************** !

subroutine S2o4_dispDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Sachin Pandey
  ! Date: 04/24/17
  !

  implicit none

  class(reaction_sandbox_s2o4_disp_type) :: this

! 12. Add code to deallocate contents of the s2o4_disp object

end subroutine S2o4_dispDestroy

end module Reaction_Sandbox_S2o4_disp_class
