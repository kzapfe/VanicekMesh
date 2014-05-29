module qd
 use precision
 use parametros
 use variables
 use funcoes
 use myfft
 implicit none

 public :: eH
 contains

! aplicar e^H a psit
  subroutine eH()
  

!  complex(adequate),dimension(:),allocatable,intent(inout) :: psio
  integer :: i
   
  !e^V psi
  forall (i=1:nx) psit(i)=eV( qs(i) )*psit(i)
!  !FT [e^V psi]
   call fft_1d(psit,nx,1)
!  !e^K (FT [e^V psi])
  forall (i=1:nx) psit(i)=eK( ps(i) )*psit(i)
!  !invFT [ e^K (FT [e^V psi]) ]
  call fft_1d(psit,nx,-1)
  !e^V invFT [ e^K (FT [e^V psi]) ]
  forall (i=1:nx) psit(i)=eV( qs(i) )*psit(i)
   

  end subroutine eH
!____________________________________________________________________________________________________
! psi_t
  subroutine evolucion()
  integer ::j 
  real(adequate)::hbar=0.1
  tau = IU*dt/(P*hbar)
  
  !do j=1,P
      call eH() 
  ! end do
  !! print *, psit 
 
  end subroutine evolucion

end module qd
