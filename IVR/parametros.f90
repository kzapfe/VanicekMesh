MODULE parametros
  USE precision
  IMPLICIT NONE
  

  COMPLEX(adequate), PARAMETER :: IU = (0.0_adequate, 1.0_adequate) ! Imaginary unit.
  REAL(adequate), parameter:: pi =3.14159265
!_____________________________________________________________________________________________________________
!CONTROL VARIABLES USED IN THE NAMELIST
!_____________________________________________________________________________________________________________
!   namelist: in_gen 
    INTEGER:: P            ! numero de splittings
    Real(adequate) :: dt            ! Number of dimensions of physical space
    REAL(adequate):: m     ! mass
    COMPLEX(adequate),dimension(:),allocatable :: psio, psit !Las funciones de onda 
    REAL(adequate),dimension(:),allocatable :: ps, qs !Momento y posicion, espero.
!   namelist: in_mc
  
! Discretizacion de Psi
    integer:: nx           ! numero de pontos da discretizacao de Psi
    real(adequate):: qmin         ! Minimo valor de qmin
    real(adequate):: dx           ! Valor del paso de la discretizacion de Psi
!   Filename
    character (LEN = 255) :: filename ! Name of the input file
END MODULE parametros
