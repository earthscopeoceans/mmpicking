real*8 function timediff(t1,t2)

! returns number of seconds t2-t1 between t1 < t2
! where t1 and t2 are integer arrays with year, julian day,
! hour, minute, seconds, milliseconds (as in SAC)

integer :: t1(6),t2(6)
integer :: m(12),kdays,n,m1,d1,m2,d2

data m/31,28,31,30,31,30,31,31,30,31,30,31/

call julian2date(t1(2),t1(1),m1,d1)
call julian2date(t2(2),t2(1),m2,d2)
kdays=march1days(t2(1),m2,d2)-march1days(t1(1),m1,d1)

timediff=kdays*86400.0+(t2(3)-t1(3))*3600.+(t2(4)-t1(4))*60.+  &
  t2(5)-t1(5)+(t2(6)-t1(6))*0.001

return
end

subroutine jul(year,month,day,jday)

! gives Julian day if called with year, month, day

implicit none

integer :: mday(12),kday(12),i,k,m,n,jday,leap,year,month,day
data mday/31,28,31,30,31,30,31,31,30,31,30,31/

if(month<1.or.month>12.or.day<1.or.day>31) then
  print *,'jul give date error ',year,month,day
  stop
endif  

leap=0
mday(2)=28
if(mod(year,4).eq.0) leap=1
if(mod(year,100).eq.0) leap=0
if(mod(year,400).eq.0) leap=1
if(leap.eq.1) mday(2)=29
kday(1)=0
do i=2,12
  kday(i)=kday(i-1)+mday(i-1)
enddo

jday=kday(month)+day

return
end

subroutine julian2date(jday,year,month,day)

! Converts Julian day jday and year into month,day

implicit none
character(len=12) :: arg

integer :: mday(12),kday(12),i,k,m,n,day,jday,month,year,leap
data mday/31,28,31,30,31,30,31,31,30,31,30,31/

if(jday>366) stop 'julian2date error in jday'
leap=0
mday(2)=28
if(mod(year,4).eq.0.and.(mod(year,100).ne.0.or.mod(year,400).eq.0))  &
  mday(2)=29
m=1
day=jday
do while(m.lt.12.and.day-mday(m)>0)
  day=day-mday(m)
  m=m+1
enddo
month=m

return
end


integer function march1days(y,m,d)

! returns integer number of days since March 1, 1900.
! (March 1, 1900 gives zero)
! method: http://alcor.concordia.ca/~gpkatch/gdate-method.html

integer :: y,m,d,mm,yy,dsincemarch

if(y<1900 .or. y>2099) then
  write(6,*) 'march1days is called with ',y,m,d
  stop 'Cannot handle year<1900 or >2099'
endif
mm=mod(m+9,12)                  ! set mm: March=0, April=1,...,Feb=11
dsincemarch=(mm*306+5)/10       ! days since March 1
yy=y-1900-mm/10                 ! yy equals year-1900 unless Jan, Feb 
march1days=365*yy+yy/4-yy/100+dsincemarch+(d-1)
if( (y.eq.2000.and.mm.gt.2) .or. (y.gt.2000) ) march1days=march1days+1

return
end

subroutine getangle(lat1,lon1,lat2,lon2,lat3,lon3,a,d12,d23)

! input: three lat/lon locations (in deg)
! finds the angle (in a flat-Earth approximation) between the
! two legs 1-2 and 2-3 in degrees
! returns distances d12 and d23 in km

implicit none
real*4, intent(in) :: lat1,lon1,lat2,lon2,lat3,lon3
real*4, intent(out) :: a,d12,d23
real*4 :: baz,d13,s,cosA,rad=0.017453292,d2km=111.194
integer :: n
logical :: db=.false.

if(db) write(13,'(a,6f10.4)') 'getangle: ',lat1,lon1,lat2,lon2,lat3,lon3
call del(lat1,lon1,lat2,lon2,d12,baz)
call del(lat1,lon1,lat3,lon3,d13,baz)
call del(lat2,lon2,lat3,lon3,d23,baz)
d12=d12*d2km
d13=d13*d2km
d23=d23*d2km
if(d12*d23>0.) then
  cosA=(d12**2+d23**2-d13**2)/(2*d12*d23)  ! cosine rule
  if(db) write(13,'(a,3f10.2,f10.5)') 'dist,cosA:',d12,d13,d23,cosA
  cosA=sign(min(abs(cosA),1.0),cosA)       ! beware of roundoff>1
  a=acos(cosA)/rad
else
  a=0.
  if(db) write(13,'(a,3f10.2,a)') 'dist,cosA:',d12,d13,d23,'    ------'
endif

return
end

subroutine del(stla,stlo,evla,evlo,dist,baz)

! input: latitudes, longitudes in degree
! dist is epicentral distance (deg) between the two points
! baz is the angle between N and the direction of (evla,evlo) as
!   seen from (stla,stlo). It is between -180 and +180, positive over East.

implicit none
real*4, intent(in) :: stla,stlo,evla,evlo
real*4, intent(out) :: dist,baz
real*8 :: slat,slon,elat,elon
real*8 :: a,scolat,ecolat,sc1,sc2,sc3,sc4,ec1,ec2,ec3,ec4,ae,be,codelb
real*8 :: halfpi=1.570796326795,pi=3.14159265359,rad=0.01745329252
real*8 :: colat
real*8 :: azi1,azi2,as,bs
logical :: db=.false.

COLAT(A)=halfpi-ATAN(.993277*TAN(A))

if(abs(stla-evla)+abs(stlo-evlo)<1.0e-5) then
  dist=0.
  baz=0.
  return
endif  

slat=stla*rad
slon=stlo*rad
elat=evla*rad
elon=evlo*rad

SCOLAT=COLAT(SLAT)
ECOLAT=COLAT(ELAT)
if(abs(slat)>halfpi) scolat=halfpi-slat
if(abs(elat)>halfpi) ecolat=halfpi-elat
SC1=SIN(SCOLAT)
SC2=COS(SCOLAT)
SC3=SIN(SLON)
SC4=COS(SLON)
EC1=SIN(ECOLAT)
EC2=COS(ECOLAT)
EC3=SIN(ELON)
EC4=COS(ELON)
AE=EC1*EC4
BE=EC1*EC3

codelb=sc1*(ae*sc4+be*sc3)+sc2*ec2
if(codelb.ge.1.0) then
  dist=0.
else if(codelb.le.-1.0) then
  dist=pi
else
  dist=acos(codelb)
endif

!as=sc1*sc4
!bs=sc1*sc3
!azi1=(as-ec3)**2+(bs+ec4)**2+sc2*sc2-2.
!azi2=(as-ec2*ec4)**2+(bs-ec2*ec3)**2+(sc2+ec1)**2-2.

azi1=(ae-sc3)**2+(be+sc4)**2+ec2*ec2-2.
azi2=(ae-sc2*sc4)**2+(be-sc2*sc3)**2+(ec2+sc1)**2-2.
if(abs(azi2).lt.1.0d-10) then
  baz=pi-sign(halfpi,azi1)
else
  baz=atan2(azi1,azi2)
endif

dist=dist/rad           ! epicentral distance in degrees
baz=baz/rad             ! back azimuth from station to event

if(db) write(6,'(a,3f10.5,f10.1)') 'lat1,lat2,cod,d=',scolat/rad,ecolat/rad, &
  codelb,dist

return
end

subroutine date2epoch(year,jday,hour,minut,nsec,epoch)

! Compute epoch (or atomic) time from date using Julian days
! The  epoch iwe use is the number of seconds that have elapsed since 
! January 1, 2018 at midnight UTC time but ignoring any leap
! seconds introduced after that.
! Incorrect before 1900 and after 2100 since it counts the century
! year also as leap year which is only true for 2000.
! tested against: https://www.epochconverter.com/

implicit none
integer, intent(in) :: year,jday,hour,minut,nsec
integer, intent(out) :: epoch
integer :: leapyears
integer*8 :: epoch8,epoch0=1514764800   ! epoch0=2018,1,0,0,0 not used
logical :: db=.false.

if(year<1970.or.year>2037) then
  print *,'date2epoch called with:',year,jday,hour,minut,nsec
  stop 'invalid year in date2epoch'
endif  

! leapyears is the number of leap years until year-1
leapyears=(year-1969)/4
epoch8 = nsec + minut*60 + hour*3600 + (jday-1)*86400 + &
  (year-1970)*31536000 + leapyears*86400

epoch=epoch8

if(db) write(13,*) 'd2ep:',leapyears,year,jday,hour,minut,nsec,epoch
return
end

subroutine epoch2date(epoch,year,jday,hour,minut,nsec)

! converts Unix epoch (TIA, seconds since Jan 1, 1970 00:00:00) to date
! (using Julian days)
! note that epoch is atomic time, i.e. UTC without any leap seconds
! tested against: https://www.epochconverter.com/

implicit none

integer, intent(in) :: epoch
integer, intent(out) :: year,jday,hour,minut,nsec
integer*8 :: epoch8,epoch0=1514764800,ns,nsl
integer :: leapyrs, ny

epoch8=epoch
! uncomment if using 2018 instead of 1970
! get Unix epoch (since 1970 instead of 2018):
! epoch8=epoch+epoch0

! subtract contribution from all years before year of current epoch
ny=1970
ns=epoch8
do while(ns>0)
  nsl=ns
  ny=ny+1
  ns=ns-31536000
  if(4*((ny-1)/4).eq.ny-1) ns=ns-86400        ! if ny is a leap year
enddo
year=ny-1
nsec=nsl
!print *,'year,nsl,nsec=',year,nsl,nsec

jday=nsec/86400+1
nsec=nsec-(jday-1)*86400
hour=nsec/3600
nsec=nsec-hour*3600
minut=nsec/60
nsec=nsec-minut*60
if(4*(year/4).eq.year) then
  if(jday>366) then
    year=year+1
    jday=1
  endif
else
  if(jday>365) then
    year=year+1
    jday=1
  endif
endif  

return
end

subroutine addt(t,tlag)
 
! t contains year,day,hour,minute,second,millisec (conform sac)
! tlag is in seconds
! tlag is added to t. 

implicit none

real*4, intent(in) :: tlag
integer, intent(inout) :: t(6)

integer :: nd
real*8 :: tsec

tsec=(t(2)-1)*86400.d0+t(3)*3600.d0+t(4)*60.d0+t(5)+t(6)*1.d-3
tsec=tsec+tlag
if(tsec.lt.0.d0) then
  t(1)=t(1)-1
  nd=365
  if(mod(t(1),4).eq.0) nd=366
  tsec=tsec+nd*86400.d0
endif
t(6)=(tsec-dint(tsec))*1000
tsec=tsec-t(6)*.001d0
t(5)=tsec-60d0*dint(tsec/60d0)
tsec=tsec-t(5)
t(4)=tsec-3600d0*dint(tsec/3600d0)
t(4)=t(4)/60.+.001
tsec=tsec-60d0*t(4)
t(3)=tsec-86400d0*dint(tsec/86400d0)
t(3)=t(3)/3600.+.001
tsec=tsec-3600.d0*t(3)
t(2)=tsec/86400d0+.001
if(t(2).lt.0) then
  t(1)=t(1)-1
  nd=365
  if(mod(t(1),4).eq.0) nd=366
  t(2)=t(2)+nd+1
  return
endif
nd=365
if(4*(t(1)/4).eq.t(1)) nd=366
if(t(2).gt.nd-1) then
  t(1)=t(1)+1
  t(2)=t(2)-nd+1
  return
endif
t(2)=t(2)+1
return
end
