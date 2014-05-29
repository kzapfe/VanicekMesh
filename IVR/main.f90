program quantum

  !! Here is only the setting of "adequate".
  use precision 
  !! aqui vienen  definidas las funciones psi y otras cosas 
  use parametros
  !! La verdad es que no sabes porque esto es un modulo
  use variables
  !! a que pedo con eso. Esto procesa el input pero no se entiende un pito.
  use input
  !! aca estan definidos los propagadores y psi inicial.
  use funcoes
  !! aqui esta la maldita transformada
  use myfft
  use qd 
  !! aqui esta el integrador de Simpson
  use num_int_
  implicit none


  integer :: i,j !! Contadores
  real(adequate)::integral
  real(adequate):: NN
  call INPUT_process_input()

  open(unit=300, file="Dinamica.dat")

 psio=psi()
 psit=psio

 write(200,*) "Esto es psio "

! do i=1,nx
!    write(300,*)qs(i), real(psio(i)), aimag(psio(i))
! enddo
 
  write(300,*)" "
  
  i=0
   
  integrand1=real(conjg(psio)*psit)
  integral=num_integration_simpson(dx,nx,integrand1)
  NN=sqrt(integral)
  write(200,*) "i'nt=",integral
  write (200,*) "Rieman suma= ", dx*sum (integrand1)
  write (200,*) "NN= ", NN
  write (100,*) i*dt,integral/(NN*NN)
  psio=psio/NN
  psit=psit/NN

  
  do i=1,P
    call evolucion()
    integrand1=real(conjg(psio)*psit)
    integral=num_integration_simpson(dx,nx,integrand1)
    NN=sqrt(integral)
     
    write (200,*) "Rieman suma= ", dx*sum (integrand1)
    write (100,*) i*dt,integral, NN
    !!print *, i*dt,real(integral),aimag(integral)

    do j=1,nx
       write(300,*)qs(j), real(psit(j)), aimag(psit(j))
    enddo
     write(300,*) " "
     write(300,*) " "
  end do



end program quantum
