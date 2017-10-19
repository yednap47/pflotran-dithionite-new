module Reaction_Sandbox_S2o4_o2_class

! Id: reaction_sandbox_s2o4_o2.F90, Tue 25 Apr 2017 09:49:22 AM MDT pandeys !
! Created by Sachin Pandey, LANL
! Description: Try to implement this reaction as zero-order wrt O2(aq)
!   - First attempt, use monod scaling factor (this was giving me problems)
!   - Second attempt, use conditional statment to calculate rate if O2(aq)>1.e-20
!------------------------------------------------------------------------------

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
    extends(reaction_sandbox_base_type) :: reaction_sandbox_s2o4_o2_type
! 3. Add variables/arrays associated with new reaction.
    character(len=MAXWORDLENGTH) :: name_spec_s2o4 ! S2O4--
    character(len=MAXWORDLENGTH) :: name_spec_o2 ! O2(aq)
    character(len=MAXWORDLENGTH) :: name_spec_so3 ! SO3--
    character(len=MAXWORDLENGTH) :: name_spec_so4 ! SO4--
    character(len=MAXWORDLENGTH) :: name_spec_h ! H+
    PetscInt :: id_spec_s2o4, id_spec_o2, id_spec_so3, id_spec_so4, id_spec_h
    PetscReal :: rate_constant, eps
  contains
    procedure, public :: ReadInput => S2o4_o2Read
    procedure, public :: Setup => S2o4_o2Setup
    procedure, public :: Evaluate => S2o4_o2React
    procedure, public :: UpdateKineticState => S2o4_o2KineticState
    procedure, public :: Destroy => S2o4_o2Destroy
  end type reaction_sandbox_s2o4_o2_type

  public :: S2o4_o2Create

contains

! ************************************************************************** !

function S2o4_o2Create()
  !
  ! Allocates s2o4_o2 reaction object.
  !
  ! Author: Sachin Pandey 
  ! Date: 09/08/16 
  !

  implicit none

  class(reaction_sandbox_s2o4_o2_type), pointer :: S2o4_o2Create

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(S2o4_o2Create)
  S2o4_o2Create%name_spec_s2o4 = ''
  S2o4_o2Create%name_spec_o2 = ''
  S2o4_o2Create%name_spec_so3 = ''
  S2o4_o2Create%name_spec_so4 = ''
  S2o4_o2Create%name_spec_h = ''
  S2o4_o2Create%id_spec_s2o4 = 0
  S2o4_o2Create%id_spec_o2 = 0
  S2o4_o2Create%id_spec_so3 = 0
  S2o4_o2Create%id_spec_so4 = 0
  S2o4_o2Create%id_spec_h = 0
  S2o4_o2Create%rate_constant = 0.d0
  S2o4_o2Create%eps = 0.d0
  nullify(S2o4_o2Create%next)

end function S2o4_o2Create

! ************************************************************************** !

subroutine S2o4_o2Read(this,input,option)
  !
  ! Reads input deck for s2o4_o2 reaction parameters (if any)
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

  class(reaction_sandbox_s2o4_o2_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,S2O4_O2')
    call StringToUpper(word)

    select case(trim(word))

      ! S2o4_o2 Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     S2O4_O2
      !       S2O4_O2_INTEGER 1
      !       S2O4_O2_INTEGER_ARRAY 2 3 4
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
                           'CHEMISTRY,REACTION_SANDBOX,S2O4_O2')
      case('EPS')
        call InputReadDouble(input,option,this%eps)
        call InputErrorMsg(input,option,'EPS', &
                       'CHEMISTRY,REACTION_SANDBOX,S2O4_O2')
! 8. Inform the user if the keyword is not recognized
      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,S2O4_O2',option)
    end select
  enddo

end subroutine S2o4_o2Read

! ************************************************************************** !

subroutine S2o4_o2Setup(this,reaction,option)
  !
  ! Sets up the s2o4_o2 reaction either with parameters either
  ! read from the input deck or hardwired.
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_s2o4_o2_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize
  this%name_spec_s2o4 = 'S2O4--'
  this%name_spec_o2 = 'O2(aq)'
  this%name_spec_so3 = 'SO3--'
  this%name_spec_so4 = 'SO4--'
  this%name_spec_h = 'H+'

  this%id_spec_s2o4 = &
    GetPrimarySpeciesIDFromName(this%name_spec_s2o4,reaction,option)
  this%id_spec_o2 = &
    GetPrimarySpeciesIDFromName(this%name_spec_o2,reaction,option)
  this%id_spec_so3 = &
    GetPrimarySpeciesIDFromName(this%name_spec_so3,reaction,option)
  this%id_spec_so4 = &
    GetPrimarySpeciesIDFromName(this%name_spec_so4,reaction,option)
  this%id_spec_h = &
    GetPrimarySpeciesIDFromName(this%name_spec_h,reaction,option)

end subroutine S2o4_o2Setup

! ************************************************************************** !

subroutine S2o4_o2React(this,Residual,Jacobian,compute_derivative, &
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

  class(reaction_sandbox_s2o4_o2_type) :: this
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

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! L_water from m^3_water

  ! mole/l-s
  ! rate = 0.0
  ! if (rt_auxvar%total(this%id_spec_o2,iphase)>1e-20) then
  !   rate = this%rate_constant*(rt_auxvar%pri_molal(this%id_spec_s2o4)+1.d-20)**1.5
  ! endif

  rate = this%rate_constant*(rt_auxvar%pri_molal(this%id_spec_s2o4))*rt_auxvar%pri_molal(this%id_spec_o2)

  if (rate < this%eps) then
    rate = 0.d0
  endif
 
! for debug
#if 0
 print *, 'k_dS2O4dt', this%k_dS2O4dt
 print *, 's2o4', rt_auxvar%total(this%id_spec_s2o4,iphase)
 print *, 'FeIII', FeIII
 print *, 'hco3', rt_auxvar%total(this%id_spec_hco3,iphase)
 print *, 'dS2O4dt', dS2O4dt
 print *, 'dHCrO4dt', dHCrO4dt
#endif

  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second  
  Residual(this%id_spec_s2o4) = Residual(this%id_spec_s2o4) - (-1.0*rate) * L_water
  Residual(this%id_spec_o2) = Residual(this%id_spec_o2) - (-1.0*rate) * L_water
  Residual(this%id_spec_so3) = Residual(this%id_spec_so3) - (1.0*rate) * L_water
  Residual(this%id_spec_so4) = Residual(this%id_spec_so4) - (1.0*rate) * L_water
  Residual(this%id_spec_h) = Residual(this%id_spec_h) - (2.0*rate) * L_water

  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

     option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
                        'due to assumptions in S2o4_o2'
     call printErrMsg(option)

  endif

end subroutine S2o4_o2React

! ************************************************************************** !

subroutine S2o4_o2KineticState(this,rt_auxvar,global_auxvar, &
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

  class(reaction_sandbox_s2o4_o2_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  ! the following arrays must be declared after reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1

end subroutine S2o4_o2KineticState

! ************************************************************************** !

subroutine S2o4_o2Destroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  implicit none

  class(reaction_sandbox_s2o4_o2_type) :: this

! 12. Add code to deallocate contents of the s2o4_o2 object

end subroutine S2o4_o2Destroy

end module Reaction_Sandbox_S2o4_o2_class
