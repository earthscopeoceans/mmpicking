program picklist

! compile: g7 picklist delaz sac bin

! command that returns a list of stations starting with the smallest
! epicentral distance and then the next closest station, etc.
! The idea is that the waveform changes little for closely spaced
! stations, making picking easier

! To create files for event/stations (including pstext files), use listsac

parameter(NSAC=8)     ! dummy dimension since data are never used
parameter(NST=500)

character*50 outfile
character*50 sacfile
character*30 date

dimension d(NSAC),dist(NST),rlat(NST),rlon(NST),list(NST)
character*8 kstnm,kvar1,kvar2
character*16 kztime
character*50 stations(NST)

jdebug=0

n=command_argument_count()    ! get number of arguments (sac files) in call
if(n.le.1) then
  print *,'Usage: picklist file1,file2,file3,...'
  print *,'e.g.: picklist *.BHZ'
  stop
endif  

! read files
k=0
kclosest=0
dclosest=99999.
do i=1,n
  call getarg(i,sacfile)
  call rsac1(sacfile,d,npts,beg,del,NSAC,nerr)
  if(nerr==-803) nerr=0           ! not a true error (sac file size)
  if(nerr /= 0) then
    print *,'SAC file read error',nerr
    print *,'in file #',i,': ',sacfile
    stop
  endif
  k=k+1
  stations(k)=trim(sacfile)
  if(k>NST) stop 'Cannot handle more than NST stations at a time'
  ! get SAC header variables, but compute gcarc in case unknown
  call getkhv('kstnm',kstnm,nerr)
  call getnhv('nzyear',nzyear,nerr)
  call getnhv('nzjday',nzjday,nerr)
  call getnhv('nzhour',nzhour,nerr)
  call getnhv('nzmin',nzmin,nerr)
  call getnhv('nzsec',nzsec,nerr)
  call getnhv('nzmsec',nzmsec,nerr)
  write(kztime,'(i2.2,a1,i2.2,a1,i2.2,f2.1)') nzhour,':',nzmin,':',  &
    nzsec,0.001*min(900,nzmsec)
  call getfhv('stla',stla,nerr)
  call getfhv('stlo',stlo,nerr)
  call getfhv('evla',evla,nerr)
  call getfhv('evlo',evlo,nerr)
  call delaz(stla,stlo,evla,evlo,DEL,AZIS,AZIE)
  dist(k)=del
  rlat(k)=stla
  rlon(k)=stlo
  if(k.eq.1) then
    slat=evla
    slon=evlo
  endif  
  if(del<dclosest) then
    kclosest=k
    dclosest=del
  endif  
  if(jdebug>0) write(13,'(i5,1x,a16,f8.2,i5,f8.2)') k,trim(sacfile),del,  &
     kclosest,dclosest
enddo  

! now walk through all stations, pick closest to the last one added to list
k=kclosest
dist(k)=-dist(k)
i=1
list(1)=kclosest
stla=rlat(kclosest)     ! save lat/lon of last station in the list
stlo=rlon(kclosest)
do while(i<n-1)
  i=i+1
  dclosest=99999.
  if(jdebug>0) write(13,*) 'Search station closest to ',trim(stations(kclosest))
  do j=1,n
    if(dist(j)<0.) cycle
    call delaz(stla,stlo,rlat(j),rlon(j),DEL,AZIS,AZIE)
    if(del<dclosest) then
      kclosest=j
      dclosest=del
      if(jdebug>0) write(13,'(2i5,1x,a16,f8.2,a2)') i,j,trim(stations(j)),del,'<-'
    else  
      if(jdebug>0) write(13,'(2i5,1x,a16,f8.2)') i,j,trim(stations(j)),del
    endif
  enddo
  list(i)=kclosest
  dist(kclosest)=-dist(kclosest)
  stla=rlat(kclosest)     ! save lat/lon of last station in the list
  stlo=rlon(kclosest)
  if(jdebug>0) write(13,'(2i5,f8.2,1x,a16)') i,kclosest,dclosest,  &
     trim(stations(kclosest))
enddo  
do i=1,n
  if(dist(i)>0.) then
    list(n)=i
    exit
  endif  
enddo  

if(jdebug>0) then
  do j=1,n
    i=list(j)
    write(13,*) i,dist(i),stations(i)
    flush(13)
  enddo  
endif

! Response has linefeeds for every 500 stations
write(6,'(500a)') (trim(stations(list(i)))//' ',i=1,n)

end

