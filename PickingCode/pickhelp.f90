program pickhelp

! g7 pickhelp ttak135 sac bin

! Use this to plot SPPIM seismograms 
! Plot seismograms in SAC format 
! Call: pickhelp file1.sac,file2.sac,...

! All files should be of the same event
! Plot also the seismograms in the order of the arguments 
! to call this program with a good order of distance and azimuth try:
! set a = `picklist 2*.sac *.SAC`
! pickhelp $a

use ttak135

implicit none
character*80 :: fname(100),dumf
character :: txt*30,clon*6,clat*5,color*11,knetwk(100)*8
character*5 :: kstnm(100)
integer, parameter :: NSAC=80000
integer :: i,i1,i2,ios,j,k,m,n,nn,nerr,npts(100)
real*4 :: y(80000,100),gcarc(100),depmax,dums(2000),x,t2x
real*4 :: lat0,lon0,beg(100),dt(100),az(100),t0(100),dx
real*4 :: stla(100),stlo(100),amp,a,r,d2r=0.01745329
real*4 :: evla,evlo,evdp,tlen,xlen,tf(100),ta(100)

! initialize
amp=1.0         ! peak ampl of plotted seismogram in cm
tlen=20.        ! duration in seconds
xlen=5.0        ! plot length in cm

n=command_argument_count()
if(n<1) stop 'Usage: pickhelp file1.sac,file2.sac,...'
print *,'pickhelp called with ',n,' files'

open(1,file='gmtpk',action='write')

do m=1,n
  call get_command_argument(m,fname(m))
  call rsac1(fname(m),y(1,m),npts(m),beg(m),dt(m),NSAC,nerr)
  if(abs(nerr).eq.803) nerr=0
  if(nerr.ne.0) then
    print *,'nerr=',nerr,' for ',trim(fname(m))
    stop
  endif  
  if(m.eq.1) then               ! get event info
    call getfhv('evla',evla,nerr)
    if(nerr.ne.0) stop 'Hypocentr absent from sac header'
    call getfhv('evlo',evlo,nerr)
    call getfhv('evdp',evdp,nerr)
    write(clat,'(f5.1)') evla
    write(clon,'(f6.1)') evlo
    clat=adjustl(clat)
    clon=adjustl(clon)
  endif
  call getkhv('kstnm',kstnm(m),nerr)
  call getfhv('gcarc',gcarc(m),nerr)
  call getfhv('az',az(m),nerr)
  call getfhv('t0',t0(m),nerr)
  call getfhv('F',tf(m),nerr)
  call getfhv('A',ta(m),nerr)
  call getkhv('knetwk',knetwk(m),nerr)
  if(t0(m)<0.) then
    t0(m)=Ptime(gcarc(m),evdp)
    if(t0(m)<0.) t0(m)=0.
    call setfhv('t0',t0(m),nerr)
    call wsac0(fname(m), x, y(1,m), nerr)
  endif  
enddo

write(1,'(a)') '#! /bin/csh'
write(1,'(a)') 'gmt set GMT_THEME modern'
write(1,'(a)') 'gmt set PROJ_LENGTH_UNIT cm'
write(1,'(a)') 'gmt set FONT_ANNOT_PRIMARY 6p'
write(1,'(a)') 'gmt set FONT_LABEL 6p'
write(1,'(a)') 'gmt set FONT_TAG 6p'         
write(1,'(a)') 'gmt set FONT_TITLE 6p'      
write(1,'(a)') 'gmt begin pplt1'
write(1,'(2a)') 'gmt subplot begin 5x3 -Fs5/4 -A+v ', &
      '-M0.4c -Scb -Srl -R-10/10/-1/1'

! write temporary scaled seismograms 
do m=1,n
  i1=(t0(m)-beg(m)-0.5*tlen)/dt(m)
  i1=max(1,i1)
  i2=i1+tlen/dt(m)
  i2=min(npts(m),i2)
  nn=i2-i1+1
  nn=min(nn,2000)
  i2=i1+nn-1
  j=0
  depmax=0.
  do i=i1,i2
    j=j+1
    dums(j)=y(i,m)
    depmax=max(abs(dums(j)),depmax)
  enddo
  dums(1:nn)=amp*dums(1:nn)/depmax
  write(dumf,'(a,i2.2)') 'dums',m
  open(2,file=dumf,action='write')
  if(nn.le.1) then
    write(6,'(2a,2f8.1,a,3i6,a)') kstnm(m),' t0,beg=', &
      t0(m),beg(m),' i1,i2,nn=',i2,i2,nn,' failed'
! else
!   write(6,'(2a,1x,i6,a)') kstnm(m),trim(dumf),nn,' samples'
  endif  
  x=-0.5*tlen
  do i=1,nn
    write(2,'(2f10.3)') x,dums(i)
    x=x+dt(m)
  enddo
  write(txt,'(a5,a,f6.1,f7.1)') kstnm(m),' @~D@~,az=',gcarc(m),az(m)
  write(1,'(3a)') '  gmt subplot set -A"',trim(txt),'"'
  color='black'
  if(knetwk(m)(1:2).eq.'MH') color='red2'
  write(1,'(4a)') '  gmt plot ',trim(dumf),' -Wthinnest,',color
  write(1,'(a)') '  gmt plot -Wthin,green << eof'
  write(1,'(a)') '0 -0.4'
  write(1,'(a)') '0 +0.4'
  if(tf(m)>-999.) then
    write(1,'(a)') '> -Wthin,purple'
    write(1,'(2f8.2)') tf(m)-t0(m),-0.4
    write(1,'(2f8.2)') tf(m)-t0(m),+0.4
  endif
  if(ta(m)>-999.) then
    write(1,'(a)') '> -Wthin,blue'
    write(1,'(2f8.2)') tf(m)-t0(m),-0.4
    write(1,'(2f8.2)') tf(m)-t0(m),+0.4
  endif
  write(1,'(a)') 'eof'
  close(2)
  if(mod(m,15).eq.0) then
    write(1,'(a)') 'gmt subplot end'
    !write(1,'(a)') 'gmt end show'
    write(1,'(a)') 'gmt end'
    write(1,'(a,i1)') 'gmt begin pplt',m/15+1
    if(m.ne.n) write(1,'(2a)') 'gmt subplot begin 5x3 -Fs5/4 -A+v ', &
      '-M0.4c -Scb -Srl -R-10/10/-1/1'
  endif  
enddo  

if(mod(n,15).ne.0) then
  write(1,'(a)') 'gmt subplot end'
  !write(1,'(a)') 'gmt end show'
  write(1,'(a)') 'gmt end'
endif
close(1)
call system('chmod a+x gmtpk')
call system('gmtpk')

end

