subroutine splineinterp(x,xa,ya,n,y2a,y)
!-------------------------------------------------------------------------------
!real*8
!	Input:
!				x--> independent variable that you want the interpolated value for
!       (i.e. altitude or energy)
!				Type: real*8
!
!				xa--> original independent variable (i.e. altitude or energy)
!				Type: real*8 1-D array
!
!				ya--> original dependent variable (i.e. density or cross-section)
!				Type: real*8 1-D array
!
!				n--> size of arrays
!				Type: integer
!
!				y2a--> second derivative of the dependent variable that you want interpolated
!       (i.e. density or cross-section)
!				Type: real*8 1-D array
!
!	Output:
!				y--> interpolated value of the dependent variable
!				Type: real*8
!
!-------------------------------------------------------------------------------
!this subroutine uses the coefficients calculated in subroutine spline. Code used
!is from Press, description from Press is listed below
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with
!xai's in order), and given the array y2a(1:n), which is the output from spline
!above, and given a value of x, this routine returns a cubic-spline interpolated
!value y
!-------------------------------------------------------------------------------------
integer n
real*8 x,y, dydx, xa(n), y2a(n), ya(n)
integer k, khi, klo
real*8 a, b, h

klo = 1
khi = n
!write(*,200) x!xa,ya,n,y2a,y
1 if(khi-klo.gt.1) then
    k = (khi+klo)/2

    if(xa(k).gt.x) then
      khi = k
    else
      klo = k
    endif
    goto 1
  endif

h = xa(khi) - xa(klo)
!print*, 'h is: ', h

if(h.eq.0.) then
  print*, 'bad xa input in splineinterp'
  stop
endif
a = (xa(khi) - x)/h
b = (x - xa(klo))/h
!print*, 'This is a and b: '
!print*, a, b

y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6

!-------------------------------------------------------------------------------------
!xa(n) contains the x-values
!ya(n) has the function values at xa(n)
!y2a(n) has the second derivatives
!klo = j and khi = j+1
!-------------------------------------------------------------------------------------

dydx = (ya(khi) - ya(klo))/h - ((3*a**2-1)*h*y2a(klo))/6 + (3*b**2-1)*h*y2a(khi)/6

return
end subroutine
