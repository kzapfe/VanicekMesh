module num_int_
  use precision
  implicit none
  
  
  public :: &
       num_integration_trapezoidal, &
       num_integration_simpson, &
       num_integration_simpson_3db8, &
       num_integration_boole

!  private ::
contains
  
  function num_integration_trapezoidal(h,npoints,funct_vals) result(output)
    ! real(adequate), external :: functiontoint
    real(adequate) :: f_0, f_1,h
    real(adequate) :: funct_vals(npoints)
    real(adequate) :: output
    integer :: i, npoints
    
    output  = 0.0
    f_0 = 0.0
    f_1 = 0.0
    
    do i = 2, npoints
       f_0 = f_0 + funct_vals(i-1)
       f_1 = f_1 + funct_vals(i)
    end do
    output = 0.5 * h * ( f_0 + f_1 )

    !! @todo to take into account an even number of points
    
  end function num_integration_trapezoidal


pure function num_integration_simpson(h,npoints,funct_vals) result(output)
    ! real(adequate), external :: functionToInt
    real(adequate),intent(in)    :: h
    complex(adequate) :: f_0, f_1,f_2
    real(adequate),intent(in)  :: funct_vals(npoints)
    real(adequate) :: output
    integer,intent(in)::  npoints
    integer :: i
    
    output  = 0.0
    f_0 = 0.0
    f_1 = 0.0
    f_2 = 0.0
    
    do i = 2, npoints-1, 2
       f_0 = f_0 + funct_vals(i-1)
       f_1 = f_1 + funct_vals(i)
       f_2 = f_2 + funct_vals(i+1)
    end do
    output = h/3.0 * (f_0 + 4.0*f_1 + f_2)

    !! To take into account an even number of points
    if( mod(npoints,2) == 0 ) then
       output = output + ( h/12.0 * ( 5.0 * funct_vals(npoints) + 8.0*funct_vals(npoints-1) - funct_vals(npoints-2) ) )
    end if
    
  end function num_integration_simpson


  function num_integration_simpson_3db8(h,npoints,funct_vals) result(output)
    ! real(adequate), external :: functionToInt
    real(adequate) :: f_0, f_1,f_2,f_3, h
    real(adequate) :: funct_vals(npoints)
    real(adequate) :: output
    integer :: i, npoints
    
    output  = 0.0
    f_0 = 0.0
    f_1 = 0.0
    f_2 = 0.0
    f_3 = 0.0

    do i = 2, npoints-2, 3
       f_0 = f_0 + funct_vals(i-1)
       f_1 = f_1 + funct_vals(i)
       f_2 = f_2 + funct_vals(i+1)
       f_3 = f_3 + funct_vals(i+2)
    end do
    output = 0.375 * h  * ( f_0 + 3.0*(f_1 + f_2) + f_3 )

    !! @todo: to take into account an even number of points
    
  end function num_integration_simpson_3db8


  function num_integration_boole(h,npoints,funct_vals) result(output)
    ! real(adequate), external :: functionToInt
    real(adequate) :: f_0, f_1,f_2,f_3,f_4, h
    real(adequate) :: funct_vals(npoints)
    real(adequate) :: output
    integer :: i, npoints
    
    output  = 0.0
    f_0 = 0.0
    f_1 = 0.0
    f_2 = 0.0
    f_3 = 0.0
    f_4 = 0.0

    do i = 2, npoints-3, 4
       f_0 = f_0 + funct_vals(i-1)
       f_1 = f_1 + funct_vals(i)
       f_2 = f_2 + funct_vals(i+1)
       f_3 = f_3 + funct_vals(i+2)
       f_4 = f_4 + funct_vals(i+3)
    end do
    output = 2.0/45.0 * h  * &
         ( (7.0 * (f_0 + f_4)) + (32.0 * (f_1 + f_3)) + (12.0 * f_2) )

    !! @todo: to take into account an even number of points
    
  end function num_integration_boole

end module num_int_
