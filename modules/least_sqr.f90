module Least_Sqr
!Purpose:
! Implement linear regression

contains

logical function test_least_squares_line() result (res)
 implicit none
 double precision, dimension(5) :: x,y
 double precision a, b, chi2, a_expect, b_expect
 x = (/1.d0, 2.d0, 3.d0, 5.d0, 7.d0/)
 y = (/2.d0, 3.d0, 3.d0, 6.d0, 8.d0/)
 call least_squares_line(5, x, y, a, b, chi2)
 a_expect = 1.025862d0
 b_expect = 0.70689655d0
 if (abs(a - a_expect) < 1.d-6 .and. abs(b - b_expect) < 1.d-6)  then
   print *, "least squares test passed"
   res = .true.
 else
   print *, "least squares test failed:"
   print *, "a = ", a, ", expected: ", a_expect
   print *, "b = ", b, ", expected: ", b_expect
   res = .false.
 endif
end function test_least_squares_line

subroutine least_squares_line(n, x, y, a, b, chi2)
 implicit none
 integer, intent(in) :: n
 double precision, intent(in), dimension(n) :: x,y
 double precision, intent(out) :: a, b, chi2
 double precision xy_av, x2_av, x_av, y_av

 xy_av = sum(x(1:n)*y(1:n)) / n
 x2_av = sum(x(1:n)*x(1:n)) / n
 x_av  = sum(x(1:n)) / n
 y_av  = sum(y(1:n)) / n

 a = (xy_av - x_av*y_av) / (x2_av - x_av*x_av)
 b = - (xy_av * x_av - y_av * x2_av) / (x2_av - x_av*x_av)

 chi2 = sum((y - a*x - b)*(y - a*x - b))

end subroutine


end module Least_Sqr
