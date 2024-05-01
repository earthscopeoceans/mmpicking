program pnext

! Compile: gfortran -o ~/bin/pnext pnext.f90

! dopick next Qnn directory listed in the file <rundopicks> (line 14)

! example of <rundopicks>:
! dopick 2018/20180907/DATA/Q00
! dopick 2018/20180908/DATA/Q00
! dopick 2018/20180908/DATA/Q01
! etc.

character :: line*50, fname*100

fname='rundopicks'
open(1,file=fname)
open(2,file='dumpnxt',action='write')

call system('grep exit '//trim(fname)//' > dum')
open(3,file='dum',action='read')
read(3,*,iostat=ios) line
if(ios.ne.0) goto 10
  
do 
  read(1,'(a)',iostat=ios) line
  if(is_iostat_end(ios)) exit
  if(line(1:1).eq.'d') then
    write(2,'(2a)') '#',trim(line)
    exit
  else
    write(2,'(a)') trim(line)
  endif
enddo
read(1,'(a)',iostat=ios) line
if(is_iostat_end(ios)) stop
if(line(1:4).ne.'exit') stop 'exit error'

10 read(1,'(a)',iostat=ios) line
if(is_iostat_end(ios)) stop
write(2,'(a)') trim(line)
write(2,'(a)') 'exit'

do
  read(1,'(a)',iostat=ios) line
  if(is_iostat_end(ios)) exit
  write(2,'(a)') trim(line)
enddo
close(1)
close(2)
close(3)

call system('mv dumpnxt '//trim(fname) )
call system('chmod a+x '//trim(fname) )
call system('./'//trim(fname) )

end
