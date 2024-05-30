program rdapf2

! compile:  g7 rdapf2 delaz timedel sac bin
! or:
! gfortran -o rdapf2 rdapf2.f90 timedel.f90 delaz.f sacio.a

! Run this program from the directory YYYYMMDD/DATA/Q* that has the APF, 
! quake and out.cfneic files

! This is a simplified version of rdapf *** for SPPIM data analysis only ***

! Program to read the APF files for a LOC pick with macro doppk in SAC
! Run this program immediately after having created an APF file with SAC,
! from the same directory as the APF and the SAC files.

! Needs to be called from the directory where the file APF resides,
! as well as where it finds the file <quake> (from qinfo).

! Arguments: phase and estimate for the pick standard errors in sec, 
! and (optional) tbias if the pick was for a later onset e.g.:
! rdapf2 P 0.3 3.5

! Since APF does not have the phase, we assume all readings in APF are for
! the same phase (P or PKP) that is given as argument. 

! Program rdapf2.f90 converts the APF file to a file that can be 
! inverted with the tomographic software suite: raydata5, voxelmatrix, 
! mpisolvetomo.

! The lines in the APF file are like:

!II.KIV.00.BHZ.M.2020.244.173227.SAC      LOC      2020244 17 33 26.83

! If multiple readings for one station are in the APF (becuase you
! changed your mind and repeated a pick typing "A,L") only the last reading
! is used by rdapf2
! If the output file add_to_datafile.* exists already, new readings are
! appended to it (which might lead to duplicate picks in output if you
! were not careful and read the same station again!).

! The origin time is read from the file <quake>. 
! From this rdapf2 computes the travel time (absolute, in seconds),
! and constructs part of the data file needed for processing by raydata5.

implicit none

! weed out obvious blunders
real*4, parameter :: MAXDELAY=20.                       ! reject > 20s

integer :: jdebug=0
integer :: julyd               ! function

! char variables
character*8 :: stationcode(2000),netw(2000),comp(2000),kstnm(100)
character :: phase*8,dum*8,filename*50,cdum*50
character*200 :: line
character*1 :: yn
character*6 :: dum6

! time variables
integer :: t0(6),tp(6)          ! origin, pick time
integer :: iodate0,iotime0,iodate,iotime,year,jday,kyd,hour,minute
real*4 :: slon0,slat0,sdep0                      ! from quake file
real*4 :: tobs(2000),tsigd(2000)

! other
logical :: ruthere
integer :: allobs
integer :: i,j,k,m
integer :: idate
integer :: ievt
integer :: iobs
integer :: ios
integer :: kluster
integer :: kpole
integer :: kunit        ! 1=nm, 2=nm/s, 3=Pa
integer :: nband
integer :: nbt
integer :: nclean
integer :: nobsa
integer :: nobst
integer, parameter :: NSAC=200000
integer :: nerr
integer :: nreject
integer :: nsmp
integer :: nmm          ! number of mermaid records in out.cfneic
integer :: replcd(2000)

real*4 :: beg
real*4 :: corcoef
real*4 :: del
real*4 :: d(NSAC)
real*4 :: evla,evlo,evdp
real*4 :: gcarc
real*4 :: ptime
real*4 :: relev(2000),rdep,rlat(2000),rlon(2000)
real*4 :: dels,azis,azie
real*4 :: rms0
real*4 :: sec
real*4 :: sdep          ! hypocentre from event file is without 0
real*4 :: slat
real*4 :: slon
real*4 :: tbias         ! difference between ISC's To and pisk wave To
real*4 :: xcl
real*4 :: Mw
real*4 :: locerr(100)   ! time error due to float mislocation
real*4 :: locerr2
real*4 :: crerr,crerr2  ! estim error in crust correction at inversion time
real*4 :: snr(100)
real*4 :: rmsb(2000)
real*4 :: pickerr       ! estimated picking error

real*8 :: timediff

k=command_argument_count()
if(k<1) then
  print *,'Usage: rdapf2 phase pickerr [tbias], eg: rdapf2 P 0.4 3.1'
  stop
endif  
call getarg(1,phase)
call getarg(2,cdum)
read(cdum,*) pickerr
tbias=0.
if(k.eq.3) then
  call getarg(3,cdum)
  read(cdum,*) tbias
endif  

stationcode=''
netw=''
comp=''
nreject=0
crerr=0.4       ! estimate of crustal correction errors for inversion
crerr2=crerr*crerr

inquire(file='APF',exist=ruthere)
if(ruthere) then       ! remove empty lines from APF file
  call system('sed -i -e "/^$/ d" APF')
else
  print *,'Cannot find file APF in this directory'
  stop
endif  
open(1,file='APF',iostat=ios)
if(ios /= 0) stop 'Error opening file APF in this directory'

open(31,file='out.cfneic',action='read',iostat=ios)
read(31,*)      ! skip header
nmm=0
snr=0.          ! signal rms for Mermaids, not used for GSN
do 
  nmm=nmm+1
  if(nmm>100) stop 'Increase dimension for locerr etc'
  read(31,'(a200)',iostat=ios) line
  if(is_iostat_end(ios)) exit
  read(line(191:197),*,iostat=ios) locerr(nmm)
  if(ios.ne.0) stop 'Error reading locerr in out.cfneic file'
  read(line(70:74),'(a)',iostat=ios) kstnm(nmm)
  if(ios.ne.0) stop 'Error reading kstnm in out.cfneic file'
  if(jdebug>0) write(13,*) ios,nmm,kstnm(nmm),locerr(nmm)
  read(line(127:131),*,iostat=ios) snr(nmm)
  if(ios.ne.0) stop 'Error reading snr in out.cfneic file'
  snr(nmm)=0.01*snr(nmm)
  if(jdebug>0) write(13,*) 'locerr,snr:',nmm,kstnm(nmm),locerr(nmm), &
    snr(nmm)
enddo
nmm=nmm-1
close(31)


inquire(file='out.rdapf2.'//trim(phase),exist=ruthere)
if(ruthere) then
  write(6,*) 'Appending to existing file out.rdapf2.'//trim(phase)
  open(4,file='out.rdapf2.'//trim(phase),access='append')
  write(4,'(//,a,/)') 'Appended:'
else
  open(4,file='out.rdapf2.'//trim(phase),action='write')
  write(4,*) 'Delays are w.r.t AK135 for ',phase
endif  

open(30,file='quake',iostat=ios)
if(ios.ne.0) then
  write(4,*) 'Error opening file quake'
  stop 'Error opening file quake'
endif  

! read event information from file 'quake' (variables end in 0)
read(30,*) t0(1),t0(2)
read(30,*) t0(3),t0(4),t0(5),t0(6)
read(30,*) slat0,slon0
read(30,*,iostat=ios) sdep0,Mw,ievt
if(jdebug>0) write(13,'(a,2f9.2,f6.0)') 'slat0,slon0,sdep0=',slat0,  &
  slon0,sdep0
close(30)
write(4,*) 'Event #',ievt,t0(1),t0(2),t0(3),t0(4),t0(5)
write(4,fmt='(a,f8.3,f9.3,f8.3)') 'hypocentre:',slat0,slon0,sdep0
write(4,'(/,2a)') 'station      lat     lon    elev   gcarc',  &
     '    tobs    tsig'

! create output file  (append if it does already exist)
write(filename,'(2a,i5.5)') 'add_to_datafile.',trim(phase),ievt
open(2,file=filename,access='append',iostat=ios)

! set defaults
nobst=1
nobsa=0
kluster=1
replcd=0
kpole=0         ! we are limited to delta<180 degrees!!
! the following are typical for onset times:
nband=0         ! ray theory
rmsb=0.
rms0=0.         ! noise rms not used
corcoef=1.      ! n/a for picks
nbt=0           ! frequency band 0 (onset time)
xcl=0.          ! no correlation length needed


! read APF
iobs=0     ! counts number of observed data 
if(jdebug>0) write(13,*) 'Reading APF'
do
  read(1,'(a)',iostat=ios) line
  if(is_iostat_end(ios)) exit
  read(line,*) filename,dum,kyd,hour,minute,sec
  filename=adjustl(filename)
  if(jdebug>0) write(13,*) trim(filename),trim(dum),kyd,hour,minute,sec
  inquire(file=filename,exist=ruthere)
  if(.not. ruthere) then
    write(6,*) 'WARNING!! File absent: ',trim(filename)
    write(4,'(2a)') 'WARNING!! File absent: ',trim(filename)
    cycle
  endif  
  if(jdebug>0) write(13,'(3i6,2a)') iobs,k,m,' line: ',trim(line)
  year=kyd/1000
  jday=kyd-1000*year
  iotime0=hour*10000+minute*100+sec
  if(jdebug>0) write(13,*) hour,minute,sec,iotime0
  iobs=iobs+1
  if(iobs>2000) stop 'Increase array dimensions, iobs>2000'
  if(ios.ne.0) then
    print *,'Error for iobs=',iobs,' in APF:',ios
    print *,adjustl(filename),year,jday,hour,minute,sec
    stop
  endif  
  if(jdebug>0) write(13,*) 'iobs=',iobs,' adjusted file: ',trim(filename)
  tp(1)=year
  tp(2)=jday
  tp(3)=hour
  tp(4)=minute
  tp(5)=sec
  tp(6)=1000.*(sec-tp(5))       ! millisecs

  ! get station information from the sac file
  call rsac1(trim(filename),d,nsmp,beg,del,NSAC,nerr)
  if(jdebug>0) write(13,*) 'nsmp,beg,del,nerr=',nsmp,beg,del,nerr
  ! get receiver coordinates 
  call getfhv('STLA',rlat(iobs),nerr)
  call getfhv('STLO',rlon(iobs),nerr)
  call getfhv('STDP',rdep,nerr)         ! instrument depth
  if(nerr.ne.0) rdep=0.
  call getfhv('STEL',relev(iobs),nerr)        ! station elevation
  if(nerr.ne.0) relev(iobs)=0.
  relev(iobs)=0.001*(relev(iobs)-rdep)              ! instrument elevation
  ! call getfhv('GCARC',gcarc,nerr)
  ! if(nerr.ne.0) gcarc=0.
  if(jdebug>0) write(13,*) 'Receiver location ',rlat(iobs),rlon(iobs),relev(iobs)
  call delaz(rlat(iobs),rlon(iobs),evla,evlo,gcarc,azis,azie)
  call getkhv('KNETWK',netw(iobs),nerr)
  if(jdebug>0) write(13,*) 'netw=',trim(netw(iobs)),' nerr=',nerr
  if(nerr.ne.0) netw(iobs)='XX'
  netw(iobs)=adjustl(netw(iobs))
  netw(iobs)(8:8)=' '
  call getkhv('KCMPNM',comp(iobs),nerr)
  if(jdebug>0) write(13,*) 'comp=',trim(comp(iobs)),' nerr=',nerr
  if(nerr.ne.0) comp(iobs)='XXX'
  comp(iobs)=adjustl(comp(iobs))
  comp(iobs)(8:8)=' '
  call getkhv('KSTNM',stationcode(iobs),nerr)
  if(jdebug>0) write(13,*) 'stationcode=',trim(stationcode(iobs)),  &
    ' nerr=',nerr
  if(nerr.ne.0) then
    write(6,*) 'WARNING! Station code absent in ',trim(filename),  &
      ' and set to XXXX'
    stationcode(iobs)='XXXX'
  endif
  stationcode(iobs)=adjustl(stationcode(iobs))  ! remove any blanks at start
  stationcode(iobs)(8:8)=' '
  call getfhv('T0',ptime,nerr)      ! get predicted P or PKP arrival time
  if(nerr.ne.0) ptime=9999.         ! flag if unknown

  locerr2=0.
  ! if mermaid there is a location error
  do i=1,nmm
    if(jdebug>0) write(13,*) i,kstnm(i).eq.stationcode(iobs),kstnm(i)
    if(kstnm(i).eq.stationcode(iobs)) then
      locerr2=locerr(i)**2
      rmsb(iobs)=snr(i)
      exit
    endif  
  enddo  
  if(jdebug>0) write(13,*) 'locerr,rmsb=',sqrt(locerr2),rmsb(iobs)
  tsigd(iobs)=sqrt(pickerr**2+locerr2+crerr2)

  ! compute observed travel time
  tobs(iobs)=timediff(t0,tp)-tbias
  if(tobs(iobs)<0.) tobs(iobs)=tobs(iobs)+86400.
  if(jdebug>0) write(13,*) hour,minute,sec,' tobs=',tobs(iobs)

  write(4,'(a8,2f8.2,f8.3,4f8.2)') stationcode(iobs),rlat(iobs),  &
     rlon(iobs),relev(iobs),gcarc,tobs(iobs),tsigd(iobs)

enddo  

allobs=iobs
nclean=0

do iobs=1,allobs
  k=0                   ! flags repeat observation for same station
  do i=iobs+1,allobs
    if(stationcode(iobs).eq.stationcode(i).and.netw(iobs).eq.netw(i)  &
      .and.comp(iobs).eq.comp(i)) then
      k=1
      replcd(iobs)=1
      if(jdebug>0) write(13,*) tobs(iobs),trim(stationcode(iobs)),  &
        netw(iobs),comp(iobs),' sets k=1 for i,iobs=',i,iobs
      write(4,'(a,f9.3,3a,f9.3)') 'tobs=',tobs(i),' in ',  &
        trim(stationcode(i)),' replaces tobs=',tobs(iobs)
    endif
  enddo
  if(jdebug>0) write(13,*) 'k=',k,' for iobs=',iobs,' in ', &
    trim(stationcode(iobs))
  if(k.eq.1) cycle      ! skip this pick since there is a later one
  nclean=nclean+1
  if(jdebug>0) then
    write(13,*) 'iobs=',iobs,' writing:'
    write(13,20) kyd,iotime0,ievt,kluster,stationcode(iobs),  &
    netw(iobs),comp(iobs),slat0,slon0,sdep0,rlat(iobs),rlon(iobs),  &
    relev(iobs),nobst,nobsa,kpole,phase,Mw
  endif
  write(2,20) kyd,iotime0,ievt,kluster,stationcode(iobs),netw(iobs),  &
    comp(iobs),slat0,slon0,sdep0,rlat(iobs),rlon(iobs),relev(iobs),  &
    nobst,nobsa,kpole,trim(phase),Mw
  kunit=2
  if(rmsb(iobs)>0.) kunit=3
  write(2,'(i2,2f10.3)') kunit,rms0,rmsb(iobs)
  write(2,*) nobst          ! nobst==1: only time observation
  write(2,'(3f9.2,i3,f9.1)') tobs(iobs),tsigd(iobs),corcoef,nbt,xcl
  write(2,*) nobsa          ! nobsa==0: no amplitude observation
enddo
20 format(3x,i7.7,2x,i6.6,i8,i4,1x,a8,1x,a8,1x,a3,2f9.3, &
  f7.1,2f9.3,f7.3,3i4,1x,a,f6.1)

write(4,'(a,i5)') 'Rejected by MAXDELAY threshold: ',nreject
write(4,'(a,i5)') 'Accepted: ',nclean

end

integer function julyd(iyr,imon,iday)

! returns julian data in format YYYYDDD
! input year,month,day. 
! See code for how I deal with truncated iyr. This won't work
! after 2050....

dimension mday(12)
data mday/0,31,59,90,120,151,181,212,243,273,304,334/

jj=iyr
! is iyr truncated to two digits? Assume 1950<jj<2050
if(iyr.lt.50) jj=jj+2000
if(iyr.lt.100) jj=jj+1900

! leap year?
leap=0
if(mod(jj,4).eq.0) leap=1
if(mod(jj,100).eq.0) leap=0
if(mod(jj,400).eq.0) leap=1

julyd=jj*1000+mday(imon)+iday
if(imon.gt.2) julyd=julyd+leap

return
end function julyd

