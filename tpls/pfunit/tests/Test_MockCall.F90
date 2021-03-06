#include "reflection.h"

module Test_MockCall_mod
   use TestSuite_mod
   use MockCall_mod
   implicit none
   private

   public :: suite

contains

#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))

   function suite()
      use TestSuite_mod, only: newTestSuite, TestSuite
      use TestMethod_mod, only: newTestMethod
      type (TestSuite) :: suite

      suite = newTestSuite('Test_MockCall')

      ADD(testExpectOneIntegerArgument)
      ADD(testFailExpectOneIntegerArgument)

   end function suite

   subroutine testExpectOneIntegerArgument
      use Assert_mod
      type (MockCall) :: mCall
      class (*), pointer :: q
      mCall = newMockCall('methodName')
      call mCall%expect(1)
      q => mCall%getExpectedValue()
      select type (p => q)
      type is (integer)
         call assertEqual(p, 1)
      end select

   end subroutine testExpectOneIntegerArgument

   subroutine testFailExpectOneIntegerArgument
      use Assert_mod
      type (MockCall) :: mCall
      class (*), pointer :: q

      integer, target :: one

      one = 1

      mCall = newMockCall('methodName')
      call mCall%expect(one)
      q => mCall%getExpectedValue()
      select type (p => q)
      type is (integer)
         call assertEqual(p, 2)
      end select
      call assertExceptionRaised()
   end subroutine testFailExpectOneIntegerArgument

end module Test_MockCall_mod
