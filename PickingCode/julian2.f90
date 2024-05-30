program julian

! compile: g7 julian bin

! Screen command returns Julian day if called with year, month, day
! or returns month, day if called with year, julianday

implicit none
character(len=12) :: arg

integer :: mday(12),kday(12),i,k,m,n,jday,jear,leap
data mday/31,28,31,30,31,30,31,31,30,31,30,31/


k=command_argument_count()

if(k<2 .or. k>4) then
  print *,'Julian can be called in two ways:'
  print *,'julian year month day: gives Julian day'
  print *,'julian year julianday: gives month day'
  stop
endif  

call get_command_argument(1,arg)
read(arg,*) jear
leap=0
mday(2)=28
if(mod(jear,4).eq.0.and.(mod(jear,100).ne.0.or.mod(jear,400).eq.0)) leap=1
if(leap.eq.1) mday(2)=29
kday(1)=0
do i=2,12
  kday(i)=kday(i-1)+mday(i-1)
enddo


if(k.eq.2) then
  call get_command_argument(2,arg)
  read(arg,*) jday
  m=1
  do while(m.lt.12.and.jday-mday(m)>0)
    jday=jday-mday(m)
    m=m+1
  enddo
  write(6,'(2i3)') m,jday
else
  call get_command_argument(2,arg)
  read(arg,*) m
  call get_command_argument(3,arg)
  read(arg,*) n
  jday=kday(m)+n
  write(6,'(i3)')jday
endif

end
