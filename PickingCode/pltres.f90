program pltres

! Compile: gfortran -o pltres pltres.f90 delaz.f90 mod_ttak135.o
! or on my Mac:
! g7 pltres delaz ttak135 bin

! make GMT plottable file of residuals in delta/azimuth plane
! from the add_to file

use ttak135

character :: fname*50
character*8 :: stationcode,netw,comp,kstnm,phase

n=command_argument_count()
if(n<1) stop 'Usage: pltres add_to*'
call get_command_argument(1,fname)
open(1,file='res.xy',action='write')
open(2,file=fname,action='read')
open(3,file='gmtres',action='write')
open(4,file='res.txt',action='write')

do
  read(2,*,iostat=ios) kyd,iotime,ievt,kluster,stationcode,netw,  &
    comp,evla,evlo,evdp,stla,stlo,stel
  if(is_iostat_end(ios)) exit
  read(2,*) kunit
  read(2,*) nobst          ! nobst==1: only time observation
  read(2,*) tobs
  read(2,*) nobsa          ! nobsa==0: no amplitude observation
  call delaz(stla,stlo,evla,evlo,DEL,AZIS,AZIE)
  tpred=Ptime(del,evdp)
  tcor=stel/5.0
  if(stel<-0.1) tcor=stel/1.5
  tpred=tpred+tcor
  write(1,'(3f8.2)') azie,del,tobs-tpred
  write(4,'(2f8.2,1x,a)') azie,del,trim(stationcode)
enddo

write(3,'(a)') '#! /bin/csh'
write(3,'(a)') 'gmt begin res'
write(3,'(a)') ''
write(3,'(a)') '  gmt set GMT_THEME modern'
write(3,'(a)') '  gmt set PROJ_LENGTH_UNIT cm'
write(3,'(a)') '  gmt makecpt -Cpolar -D -T-2/2/1 -Z'
write(3,'(a)') '  gmt plot res.xy -R0/360/0/95 -JP14+a -Ba30fg -Sc0.25 -C -Wthinnest'
write(3,'(a)') '  gmt text res.txt -N -B -F+f4p+jBL -D0.1'
write(3,'(a,i7,1x,i6.6,a)') '  echo 42 115 ',kyd,iotime,' | gmt text -N'
write(3,'(a)') "  gmt colorbar -Dx0/-0.6+w4/0.5+h -Ba -Bx+l'sec'"
write(3,'(a)') ''
write(3,'(a)') 'gmt end show'
close(3)

call system('chmod a+x ./gmtres')
call system('./gmtres')

end
