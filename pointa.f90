program pointa

! Compile: gfortran -o pointa pointa.f90 timedel.f90 mod_ttak135.o -g -w /usr/local/lib/libsacio_bsd.a

! or on my Mac: g7 pointa timedel ttak135 sac bin

! using the APF file of onset picks, add A pointers in the SAC headers
! if an F pointer is detected, write xy file with A-F

use ttak135

implicit none
integer, parameter :: NSAC=100000
character :: fname(100)*80,loc*3,kstnm(100)*8
integer :: nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec
integer :: pkyear,pkjday,pkhour,pkmin,pksec,pkmsec
integer :: i,j,k,ios,nerr,yearday,npts
integer :: t0(6),tp(6)
real*8 :: timediff
real*4 :: ta(100),y(NSAC),beg,dt,sec,tf(100)
real*4 :: x,gcarc(100),tpred,evdp(100)


open(1,file='APF',action='read')
open(2,file='/Users/auguste/data/Mermaids/SPPIM/af.xy',access='append')
open(3,file='/Users/auguste/data/Mermaids/SPPIM/out.pointa',access='append')
open(4,file='/Users/auguste/data/Mermaids/SPPIM/ddf.xy',access='append')
!write(3,'(a)') 'kstnm      gcarc      ta      tf     dta     dtf  filename'
k=0
tf=-12345.
ta=-12345.
do 
  k=k+1
  if(k>100) stop 'k>100'
  read(1,*,iostat=ios) fname(k),loc,yearday,tp(3),tp(4),sec
  if(is_iostat_end(ios)) exit
  tp(5)=sec
  tp(6)=1000*(sec-tp(5))
  tp(1)=yearday/1000
  tp(2)=yearday-1000*tp(1)
  call rsac1(fname(k),y,npts,beg,dt,NSAC,nerr)
  if(npts.ge.NSAC) stop 'npts exceeds NSAC'
  if(nerr.ne.0) then
    print *,'rsac1 read error ',nerr
    cycle
  endif  
  call getnhv('nzyear',t0(1),nerr)
  call getnhv('nzjday',t0(2),nerr)
  call getnhv('nzhour',t0(3),nerr)
  call getnhv('nzmin',t0(4),nerr)
  call getnhv('nzsec',t0(5),nerr)
  call getnhv('nzmsec',t0(6),nerr)
  if(nerr.ne.0) then
    print *,'Origin time reading failed for ',trim(fname(k))
    cycle
  endif  
  call getkhv('kstnm',kstnm(k),nerr)
  call getfhv('gcarc',gcarc(k),nerr)
  call getfhv('evdp',evdp(k),nerr)
  j=k
  do i=1,k-1
    if(kstnm(i).eq.kstnm(k)) then
      j=i
      ta(k)=-1.
      exit
    endif  
  enddo
  call getfhv('F',tf(j),nerr)
  if(nerr.ne.0.or.tf(j)<0.) tf(j)=-12345.
  ta(j)=timediff(t0,tp)
  call setfhv('A',ta(j),nerr)
  call setkhv('KA','A',nerr)
  call wsac0(fname(k),x,y,nerr)
enddo  
k=k-1

do j=1,k  
  if(ta(j)<0.) cycle
  tpred=Ptime(gcarc(j),evdp(j))
  if(tf(j)>0.) then
    write(2,'(2f8.2)') ta(j)-tpred,tf(j)-tpred
    write(3,'(a8,5f8.2,2x,a)') kstnm(j),gcarc(j),ta(j),tf(j),ta(j)-tpred, &
      tf(j)-tpred,trim(fname(j))
    write(4,'(2f8.2)') evdp(j),tf(j)-ta(j)
  else  
    write(3,'(a8,3f8.2,18x,a)') kstnm(j),gcarc(j),ta(j),ta(j)-tpred, &
      trim(fname(j))
  endif
enddo

close(2)
close(3)

end
