subroutine julian(year,jday,month,cmonth,day)

! find month/day from julian day

implicit none

integer, intent(in) :: year,jday
integer, intent(out) :: month, day
integer :: mday(12),leap,m
character*3 mon(12),cmonth

data mday/31,28,31,30,31,30,31,31,30,31,30,31/
data mon/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep', &
  'Oct','Nov','Dec'/

cmonth=''
day=0
! Pope Gregory decreed that October 4, 1582, on the Julian calendar 
! was to be followed 
! by October 15, 1582, in the newly established Gregorian calendar.
if(year.lt.1582 .or. (year.eq.1582.and.(month<10.or.  &
  (month.eq.10.and.day<15)))) then
    print *,'subroutine julian is based on the Gregorian calendar'
    print *,'it and does not handle the Julian calendar < AD 1582'
  return
endif  
leap=0
day=jday
mday(2)=28
if(mod(year,4).eq.0.and.(mod(year,100).ne.0.or.mod(year,400).eq.0)) leap=1
if(leap.eq.1) mday(2)=29
m=1
do while(m.lt.12.and.day-mday(m)>0)
  day=day-mday(m)
  m=m+1
enddo
month=m
cmonth=mon(m)

return
end
