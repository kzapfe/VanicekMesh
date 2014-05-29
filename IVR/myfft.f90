!EZ Adapted version of CVS/dynamics/myfft.f90 for PIMC
!   Only for 1d FFT and FFTW3
!   Here is used Datasize instead of N. 
MODULE myfft

  USE precision
  use parametros
  use funcoes
!  USE parameters
!#ifdef FFTW3
!  USE fftw3
!#endif
!#ifdef MPI
!  USE mympi_
!#endif

  IMPLICIT NONE

!#ifdef FFTW3
  EXTERNAL dfftw_plan_dft, dfftw_plan_dft_1d, dfftw_destroy_plan, dfftw_execute 
  include '/usr/include/fftw3.f'
!#endif

  INTEGER, PARAMETER :: wp = KIND(1.0D0)

!#ifdef FFTW3
  COMPLEX(wp), DIMENSION (:), ALLOCATABLE, PRIVATE :: in, out
  INTEGER*8, SAVE :: plan_forward, plan_backward
!#endif /* FFTW3 */

CONTAINS

  SUBROUTINE initialize_fft(Datasize)
    INTEGER :: j
    INTEGER :: alloc_status = 0
    ! EZ: Definition Datasize
    INTEGER,INTENT(IN):: Datasize
!print *,alloc_status

!#ifdef FFTW3
!    IF(myid == 0) THEN
!      WRITE (*,*) 'Initializing FFT. FFT library: FFTW3.'
!      WRITE (*,*)
!    END IF
    ALLOCATE(in(0:Datasize-1), out(0:Datasize-1), stat = alloc_status)
    IF (alloc_status /= 0) STOP 'Allocation fails for FFT.'
   
   ! EZ Only for 1D FFT
    CALL dfftw_plan_dft_1d(plan_forward,Datasize,in,out,FFTW_FORWARD,FFTW_ESTIMATE) 
    CALL dfftw_plan_dft_1d(plan_backward,Datasize,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
!#endif /* FFTW3 */

    
  END SUBROUTINE initialize_fft

  SUBROUTINE finalize_fft()
    INTEGER :: j

!#ifdef FFTW3
    DEALLOCATE(in, out)
    CALL dfftw_destroy_plan(plan_forward)
    CALL dfftw_destroy_plan(plan_backward)
!#endif 


  END SUBROUTINE finalize_fft

!****************** General purpose 1-D FFT. Typically used to get spectra (of varying length) during dynamics. All allocations and initializations are done inside the subrouine. Therefore it's slower than 1-D FFT used for dynamics.
  SUBROUTINE fft_1d(fdata, Datasize, sense) 

    COMPLEX(adequate), INTENT(INOUT), DIMENSION(:) :: fdata
    INTEGER, INTENT(IN) :: Datasize
    INTEGER, INTENT(IN) :: sense
    INTEGER:: alloc_status = 0

!#ifdef FFTW3
    COMPLEX(wp), DIMENSION (:), ALLOCATABLE :: outfdata
    INTEGER*8 :: plan=3
!#endif

!-1 for forward propagation + 1 for backward
    IF (sense == 1) THEN
      fdata = CONJG(fdata)
    END IF

!#ifdef FFTW3
    ALLOCATE(outfdata(Datasize), stat = alloc_status)
    IF (alloc_status /= 0) STOP 'Allocation fails for gen purpose FFT.'
    CALL dfftw_plan_dft_1d(plan,Datasize,fdata,outfdata,FFTW_FORWARD,FFTW_ESTIMATE) 
    CALL dfftw_execute(plan)
    fdata = 1 / SQRT(REAL(Datasize)) * outfdata
    CALL dfftw_destroy_plan(plan)
    DEALLOCATE(outfdata)
!#endif

!-1 for forward propagation + 1 for backward
    IF (sense == 1) THEN
      fdata = CONJG(fdata)
    END IF

  END SUBROUTINE fft_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A PEDAL


  SUBROUTINE mifft_1d(fpsi,Datasize,sense,delta)
  
  complex(adequate),dimension(:),intent(inout) :: fpsi 
  complex(adequate),dimension(:,:),allocatable :: ft_array 
  integer,intent(in) :: Datasize,sense
  integer:: j,i
  real(adequate),intent(in):: delta
  
  ALLOCATE(ft_array(Datasize,Datasize))
  ft_array=0.0_adequate

  !!Transformando de q a p
  if (sense > 0)  then
     do j=1,Datasize
        do i=1, Datasize
         ft_array(j,i) = exp(-IU*ps(j)*qs(i))*fpsi(i)
      end do
   end do
   
   forall (j=1:Datasize) fpsi(j)=sum(ft_array(j,:))      
   fpsi = delta*fpsi/sqrt(2.0*PI)
else
   
   !!Transformando de p a q
   do j=1,Datasize
      do i=1, Datasize
         ft_array(j,i) = exp(IU*ps(i)*qs(j))*fpsi(i)
      end do
      end do
       
      forall (j=1:Datasize) fpsi(j)=sum(ft_array(j,:))      
      fpsi = delta*fpsi/sqrt(2.0*PI)
  end if

  deallocate (ft_array)
  END SUBROUTINE mifft_1d




END MODULE myfft
