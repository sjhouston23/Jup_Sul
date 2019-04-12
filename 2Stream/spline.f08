subroutine spline(x, y, n, y2)
!-------------------------------------------------------------------------------
!
!	Input:
!				x--> independent variable (i.e. altitude or energy)
!				Type: real 1-D array
!
!				y--> dependent variable (i.e. density or cross-section)
!				Type: real 1-D array
!
!				n--> size of arrays
!				Type: integer
!
!	Output:
!				y2--> second derivative of the dependent variable
!				produced for the spline interpolator
!				Type: real
!
!!! SIDENOTE:yp1 and yp2 are the first derivatives of the first and last (nth)
!!! point of the dependent variable. They can be used as an input variable, but
!!! here I have just set them to 0.0 inside this subroutine.
!-------------------------------------------------------------------------------
integer n, NMAX
real yp1, ypn, x(n), y(n), y2(n)
parameter (NMAX = 500)
!-------------------------------------------------------------------------------------
!given arrays x(1:n) and y(1:n) containing a tabulared function i.e. yi = f(xi),
!with x1 < x2 < ::: < xN, and given values yp1, and ypn for the rst derivative of the
!interpolating function at points 1, and n, respectively, this routine returns an array
!cy2(1:n) of length n which contains the second derivatives of the interpolating function
!at the tabulated points xi. If yp1 and/or ypn are equal to 1e10 30 or larger, the routine
!is signaled to set the corresponding boundary condition for a natural apline, with
!zero second derivative on that boundary
!-------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------
!parameter NMAX is the largest anticipated value of n
!-------------------------------------------------------------------------------------

integer i, k
real p, qn, sig, un, u(NMAX)
yp1=0.0
ypn=0.0
if (yp1.gt..99e30) then 		!the lower boundary condition is set either to be natural
	y2(1) = 0
	u(1) = 0
else 							!or else to have a specied rst derivative
	y2(1) = -0.5
	u(1) = (.3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
endif

!-------------------------------------------------------------------------------------
!this is the decomposition loop of the tridiagonal algorithm. y2 and u are used for
!temporary storage of the decomposed factors
!-------------------------------------------------------------------------------------

do i = 2, n-1
	sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
	p = sig*y2(i-1)+2
	y2(i) = (sig-1.)/p
	u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
	/(x(i) - x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
enddo

if (ypn.gt..99e30) then 		!the lower boundary condition is set either to be natural
	qn = 0
	un = 0
else							!or else to have a specied rest derivative
	qn = 0.5
	un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
endif

y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)

do k = n-1, 1, -1 				!this is the backsubstitution loop of the tridiagonal algorithm
	y2(k) = y2(k)*y2(k+1)+u(k)
enddo

return
end subroutine
