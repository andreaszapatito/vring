module loop 
!#include <petsc/finclude/petscdef.h>

!  use petsc
  use run_tools
  use cons_tools
  use para_tools
  use comm_tools
  use mesh_tools
  use solu_tools
  use solp_tools
  use part_tools

  use poisson
  use poisson_solve

  use momentum_coef
  use momentum_terms
  use inverse
  use divergence

  use communicate
  implicit none
!#include <mpif.h>
  contains
    subroutine iterate (var)
      implicit none

      integer               :: ierr
      integer               :: code
      integer               :: istep
      integer               :: irkstep
      integer               :: id
      integer               :: nparticles
      type(run)             :: var

      integer ic,jc,kc

      nparticles=1

      allocate(var%prt(nparticles))
      var%con%timming0= MPI_Wtime()/60.0
      do id=1,nparticles
        call init_part(id,var%prt(id),var%msh,var%par)
      enddo
! Time advancement loop start
    do istep= 0,var%par%totstp
       call communicateeq(var)
!      do istep= 1,1000
! Update timestep
        var%par%istep = istep
        var%par%nstep = istep+var%par%nstepi
        var%par%ntime = var%par%dt+var%par%ntime
      do id=1,nparticles
!        write (*,*) 'inject',id,var%par%nstep,var%prt(id)%inj
        if (mod(var%par%nstep,var%prt(id)%inj)==0) then
                if (istep==0) call inlt_part(var%par,var%prt(id),var%msh,var%su,var%com)
        endif
      enddo
! Non-solenoidal component
        var%con%timming1= MPI_Wtime()/60.0
! RK Loop start
        do irkstep = 1,3 !var%par%nrkstep
          call communicateeq(var)
          do id=1,nparticles
           ! call iter_part(irkstep,var%prt(id),var%msh,var%com,var%su,var%par) particles need further temporal refinement 
            if (irkstep==1)  call iter_part(irkstep,var%prt(id),var%msh,var%com,var%su,var%par)
          enddo
          var%par%irkstep = irkstep
          call ucoeff(var)  
          call verbose(var,"ucoeff        called")
! coesys
          call coesys(var)  
          call verbose(var,"coesys        called")

 
          call communicateeq(var)
!        call momentum
          call hdnbc(var);  call verbose(var,"hdnbc         called")
          call hdnl1(var);  call verbose(var,"hdnl1         called")
          call commq(var,var%su%dq1,3)
          call commq(var,var%su%dq1,2)
          call hdnl2(var);  call verbose(var,"hdnl2         called")
          call commq(var,var%su%dq2,3)
          call commq(var,var%su%dq2,2)
          call hdnl3(var);  call verbose(var,"hdnl3         called")
          call commq(var,var%su%dq3,3)
          call commq(var,var%su%dq3,2)

          call communicateeq(var)
! inverse AF1
! inverse AF2
! inverse AF3

          do ic=1,var%msh%ntheta
            do jc=1,var%msh%nrc
              do kc=1,var%msh%nzc
!                var%su%q1(ic,jc,kc) = 0.001*(var%msh%zc(kc))**2 + 0.001*(var%msh%rc(jc))**2
!                var%su%q2(ic,jc,kc) = 0.001*(var%msh%zc(kc))**2 + 0.001*(var%msh%rc(jc))**2
!                var%su%q3(ic,jc,kc) = 0.001*(var%msh%zc(kc))**2 + 0.001*(var%msh%rc(jc))**2
              enddo
            enddo
          enddo

          call communicateeq(var)
          call step(var);  call verbose(var,"step          called")
! divergence
          call communicateeq(var) ;  call verbose(var,"comm          called")
          call divq(var);  call verbose(var,"divq          called")
          call communicateeq(var) ;  call verbose(var,"comm          called")

! Poisson solver
          call laplace_solve(var%sp,var%com,var%msh,var%par)    ;  call verbose(var,"laplace_solve called")
          call communicateeq(var)
            
!          call commq(var,var%sp%dph,3)
!          call commq(var,var%sp%dph,2)


! poisson
          call upsol(var)


          call communicateeq(var)
          call uppres(var)
          call communicateeq(var)
   

        enddo

! RK Loop end

       var%con%timming2= MPI_Wtime()/60.0
       var%con%timming3= MPI_Wtime()/60.0
       if (mod((var%par%nstep),var%par%expstp).eq.0)   call restart_save(var);  call verbose(var,"export solve  called") 

        if (mod((var%par%nstep),var%par%anlstp).eq.0) call analysis(var);  call verbose(var,"analysis called") 
!        if (mod((var%par%nstep-1),var%par%expstp).eq.0) call ffemexport(var,1);  call verbose(var,"ffemexp  called") 
!        if (mod((var%par%nstep-1),var%par%expstp).eq.0) call ffemexport(var,2);  call verbose(var,"ffemexp  called") 
!        if (mod((var%par%nstep-1),var%par%expstp).eq.0) call ffemexport(var,3);  call verbose(var,"ffemexp  called") 
!        if (mod((var%par%nstep-1),var%par%expstp).eq.0) call ffemexport(var,4);  call verbose(var,"ffemexp  called") 
!        if (mod((var%par%nstep-1),var%par%expstp).eq.0) call ffemexport(var,5);  call verbose(var,"ffemexp  called") 
        if (var%com%ip.eq.0) then
          write (*,"(2a)") "__________________________________________"
          write (*,"(2a)") "::Step________::Time________::Residual____::Performance_::Partitioning"
          write (*,"(a,I10,3(a,F10.6),a,I10)") " I: ",var%par%nstepi, " I: ",var%par%ntimei, " M1:",var%sp%r1," TM:",var%con%timming2-var%con%timming1," TN:",var%com%np
          write (*,"(a,I10,3(a,F10.6),a,I10)") " R: ",var%par%nstep,  " R: ",var%par%ntime,  " M2:",var%sp%r2," TP:",var%con%timming3-var%con%timming2," TT:",var%com%np_a(1)
          write (*,"(a,I10,3(a,F10.6),a,I10)") " F: ",var%par%nstepf, " F: ",var%par%ntimef, " P1:",var%su%r1," TS:",var%con%timming3-var%con%timming1," TR:",var%com%np_a(2)
          write (*,"(a,I10,3(a,F10.6),a,I10)") " E: ",var%par%expstp, " Dt:",var%par%dt,     " P2:",var%su%r2," TT:",(var%con%timming3-var%con%timming0)/float(istep)," TZ:",var%com%np_a(3)
        endif

      enddo
! Time advancement loop end
     
      call mpi_barrier(var%com%comm_0, ierr)
    
  !    call petscfinalize(ierr)
  end subroutine iterate

end module loop
