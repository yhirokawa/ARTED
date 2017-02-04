!
!  Copyright 2016 ARTED developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
module timer
  implicit none

  integer, public, parameter :: TIMER_DT_EVOLVE    = 0
  integer, public, parameter :: TIMER_HPSI         = 1
  integer, public, parameter :: TIMER_PSI_RHO      = 2
  integer, public, parameter :: TIMER_HARTREE      = 3
  integer, public, parameter :: TIMER_EXC_COR      = 4
  integer, public, parameter :: TIMER_CURRENT      = 5
  integer, public, parameter :: TIMER_TOTAL_ENERGY = 6
  integer, public, parameter :: TIMER_ION_FORCE    = 7
  integer, public, parameter :: TIMER_DT_EVOLVE_AC = 8
  integer, public, parameter :: TIMER_K_SHIFT_WF   = 9
  integer, public, parameter :: TIMER_OTHER        = 10
  integer, public, parameter :: TIMER_ALLREDUCE    = 11

  integer, public, parameter :: TIMER_CG           = 12
  integer, public, parameter :: TIMER_DIAG         = 13
  integer, public, parameter :: TIMER_SP_ENERGY    = 14
  integer, public, parameter :: TIMER_GRAM_SCHMIDT = 15

  integer, public, parameter :: TIMER_HPSI_INIT    = 16
  integer, public, parameter :: TIMER_HPSI_STENCIL = 17
  integer, public, parameter :: TIMER_HPSI_PSEUDO  = 18
  integer, public, parameter :: TIMER_HPSI_UPDATE  = 19

  integer, public, parameter :: TIMER_DYNAMICS     = 20

  public timer_initialize
  public timer_reentrance_read, timer_reentrance_write
  public timer_set, timer_reset
  public timer_begin, timer_end
  public timer_thread_begin, timer_thread_end

  public timer_show_hour, timer_show_min
  public timer_write, timer_thread_write
  public timer_get, timer_thread_get
  public get_wtime

  integer, private, parameter   :: TIMER_SIZE = 30
  real(8), private, allocatable :: log_time(:)
  real(8), private, allocatable :: log_temp(:)

  real(8), private, allocatable :: log_time_t(:,:)
  real(8), private, allocatable :: log_temp_t(:,:)

private
contains
  subroutine timer_initialize
    use omp_lib, only: omp_get_max_threads
    implicit none
    allocate(log_time(0:TIMER_SIZE - 1))
    allocate(log_temp(0:TIMER_SIZE - 1))
    allocate(log_time_t(0:TIMER_SIZE - 1, 0:omp_get_max_threads()-1))
    allocate(log_temp_t(0:TIMER_SIZE - 1, 0:omp_get_max_threads()-1))
    call timer_reset
  end subroutine

  subroutine timer_set(e,t)
    implicit none
    integer,intent(in) :: e
    real(8),intent(in) :: t
    log_time(e) = t
    log_temp(e) = 0.d0
  end subroutine

  subroutine timer_reset(e)
    implicit none
    integer,intent(in),optional :: e
    if(present(e)) then
      log_time  (e)   = 0.d0
      log_temp  (e)   = 0.d0
      log_time_t(e,:) = 0.d0
      log_temp_t(e,:) = 0.d0
    else
      log_time  (:)   = 0.d0
      log_temp  (:)   = 0.d0
      log_time_t(:,:) = 0.d0
      log_time_t(:,:) = 0.d0
    end if
  end subroutine

  subroutine timer_reentrance_read(fd)
    implicit none
    integer,intent(in) :: fd
    read(fd) log_time(0:TIMER_SIZE - 1)
    read(fd) log_temp(0:TIMER_SIZE - 1)
  end subroutine

  subroutine timer_reentrance_write(fd)
    implicit none
    integer,intent(in) :: fd
    write(fd) log_time(0:TIMER_SIZE - 1)
    write(fd) log_temp(0:TIMER_SIZE - 1)
  end subroutine

  subroutine timer_begin(e)
    implicit none
    integer,intent(in) :: e
    log_temp(e) = get_wtime()
  end subroutine

  subroutine timer_end(e)
    implicit none
    integer,intent(in) :: e
    log_time(e) = log_time(e) + get_wtime() - log_temp(e)
  end subroutine

  subroutine timer_thread_begin(e)
    use omp_lib
    implicit none
    integer,intent(in) :: e
    integer :: tid
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timer_begin(e)
    end if
    log_temp_t(e,tid) = get_wtime()
  end subroutine

  subroutine timer_thread_end(e)
    use omp_lib
    implicit none
    integer,intent(in) :: e
    integer :: tid
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timer_end(e)
    end if
    log_time_t(e,tid) = log_time_t(e,tid) + get_wtime() - log_temp_t(e,tid)
  end subroutine

  subroutine timer_show_hour(str, e)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: e
    real(8) :: time,hour
    time = log_time(e)
    hour = time / 3600
    write(*,*) str,time,'sec =',hour,'hour'
  end subroutine

  subroutine timer_show_min(str, e)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: e
    real(8) :: time,mini
    time = log_time(e)
    mini = time / 60
    write(*,*) str,time,'sec =',mini,'min'
  end subroutine

  subroutine timer_write(fd,str,e)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,e
    real(8) :: time
    time = log_time(e)
    write(fd,*) str,time,'sec'
  end subroutine

  subroutine timer_thread_write(fd,str,e)
    use omp_lib
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,e
    real(8) :: time
    integer :: i
    write(fd,*) str
    do i=0,omp_get_max_threads()-1
      time = log_time_t(e,i)
      write(fd,*) 'tid =',i,': ',time,'sec'
    end do
  end subroutine

  function timer_get(e)
    implicit none
    integer,intent(in) :: e
    real(8)            :: timer_get
    timer_get = log_time(e)
  end function

  function timer_thread_get(e,tid)
    implicit none
    integer,intent(in) :: e,tid
    real(8)            :: timer_thread_get
    timer_thread_get = log_time_t(e,tid)
  end function

  function get_wtime()
    implicit none
    real(8) :: get_wtime
    real(8) :: omp_get_wtime
    get_wtime = omp_get_wtime()
  end function
end module
