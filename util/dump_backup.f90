program main
  use iso_fortran_env
  use global_variables
  use timer
  implicit none
  integer :: iter,ixy_m
  character(len=128) :: fmt

  call timer_initialize
  open(unit=INPUT_UNIT,form='unformatted')
  call prep_backup_list(.FALSE., INPUT_UNIT)

  print *, 'iter_now, entrance_iter =',iter_now,entrance_iter

  print *, '==============================='

  print *, 'NKx,y,z     =',NKx,NKy,NKz
  print *, 'NKxyz       =',NKxyz
  print *, 'NK          =',NK
  print *, 'kbTev       =',kbTev
  print *, 'NB,NBoccmax =',NB,NBoccmax
  print *, 'NLx,y,z     =',NLx,NLy,NLz
  print *, 'NL          =',NL
  print *, 'Sym         =',Sym
  print *, 'NI,NE       =',NI,NE
  print *, 'Nscf,Ncg    =',Nscf,Ncg
  print *, 'Nt,dt       =',Nt,dt

  print *, '==============================='

  if (allocated(data_vac_Ac)) then
    file_ac_vac = './Ac_Vac.out'
    print *, 'out : ', trim(file_ac_vac)
    open(900,file=file_ac_vac, position='rewind')
    do iter=0,Nt
      write(900,"(7e26.16e3)")iter*dt,data_vac_Ac(1:3,1,iter),data_vac_Ac(1:3,2,iter)
    end do
    close(900)
  end if

  if (allocated(data_local_Ac)) then
    write(file_ac_m, "('./Ac_M',I6.6,'.out')") NXY_s
    print *, 'out : ', trim(file_ac_m)
    write(fmt,"(A,I2,A)")"(",(NXY_e-NXY_s+1)*6+1,"e26.16e3)"
    open(900, file=file_ac_m, position='rewind')
    write(900,"(A,2x,I6,2x,A,2x,I6)")"# Data of macro points",NXY_s,"-",NXY_e
    do iter=0,Nt
       write(900,fmt)iter*dt,(data_local_Ac(1:3,ixy_m,iter),data_local_jm(1:3,ixy_m,iter),ixy_m = NXY_s,NXY_e)
    end do
    close(900)
  end if
end program

subroutine err_finalize(err_message)
  use communication
  implicit none
  character(*),intent(in) :: err_message
  write(*,*) err_message
  stop
end subroutine
