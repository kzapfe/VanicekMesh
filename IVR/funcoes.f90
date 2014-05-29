MODULE funcoes

  use precision
  USE parametros
  use variables

  IMPLICIT NONE

  CONTAINS
! Potential
pure complex function V(x)
 real(adequate),intent(in) ::x
 complex ::V
 !!V=0.00
 V = x**2
end function V

! Kinetic Energy part of Hamiltonian
pure  COMPLEX FUNCTION eK(xx)

  real(adequate),intent(in)     :: xx
  COMPLEX :: eK
  
  eK = exp(-tau*xx**2/CMPLX(2.0*m))

  END FUNCTION eK

! Potential part of the Hamiltonian
pure  COMPLEX FUNCTION eV(qq)

  real(adequate),INTENT(IN)     :: qq
  COMPLEX :: eV
  
  eV =EXP(-tau*V(qq)/2.0)

  END FUNCTION eV

! ESTADO INICIAL 
  COMPLEX function psi()
  
  complex,dimension(:),allocatable :: psi
  integer :: i

  ALLOCATE(psi(nx))
  
  
  FORALL (i=1:nx) 
    psi(i)=EXP(-0.5*(qs(i)-1.0)**2)
  END FORALL
  END FUNCTION psi

  
END MODULE funcoes
