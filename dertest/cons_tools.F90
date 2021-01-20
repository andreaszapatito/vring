module cons_type
  implicit none

  type :: cons
    character(len=10)     :: inp_file
    character(len=10)     :: out_file
    character(len=10)     :: rst_file
    character(len=10)     :: fld_file

    real(kind=8)          :: pi

    real(kind=8)          :: timming0

    real(kind=8)          :: timming1
    real(kind=8)          :: timming2
    real(kind=8)          :: timming3
    real(kind=8)          :: timming4
    real(kind=8)          :: timming5
    real(kind=8)          :: timming6
    real(kind=8)          :: timming7
    real(kind=8)          :: timming8


    real(kind=8)          :: timmingPO
    real(kind=8)          :: timmingMO
    real(kind=8)          :: timmingC1
    real(kind=8)          :: timmingC2

    real(kind=8)          :: rk_al
    real(kind=8)          :: rk_ga
    real(kind=8)          :: rk_ro
    


  end type cons

end module cons_type

module cons_tools
  use cons_type

  implicit none

contains

  subroutine create_cons(con)
    implicit none
    type(cons)         :: con


    con%pi=4.d0*atan(1.0)
    con%inp_file='inp.txt'  
    con%out_file='out.txt'  
    con%rst_file='rst.dat'  
    con%fld_file='fld'  
 
  end subroutine create_cons

end module cons_tools
