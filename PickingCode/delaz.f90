subroutine delaz(stla,stlo,evla,evlo,DEL,AZIS,AZIE)

! this subroutine is needed in several programs
! in contrast to azidl, this routine works with DEGREE input/output


COLAT(A)=1.57079632679d0-ATAN(.993277d0*TAN(A))

pi=3.14159265359
rad=180./pi

if(abs(stla-evla)+abs(stlo-evlo)<1.0e-5) then
  dist=0.
  baz=0.
  return
endif  

slat=stla/rad
slon=stlo/rad
elat=evla/rad
elon=evlo/rad
SCOLAT=COLAT(SLAT)
ECOLAT=COLAT(ELAT)
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
AZI1=(AE-SC3)**2+(BE+SC4)**2+EC2*EC2-2.
AZI2=(AE-SC2*SC4)**2+(BE-SC2*SC3)**2+(EC2+SC1)**2-2.
IF(abs(AZI2)<1.0e-10) then
  AZIS=3.141592653-SIGN(1.570796326,AZI1)
else
  AZIS=ATAN2(AZI1,AZI2)
endif
CODELB=SC1*(AE*SC4+BE*SC3)+SC2*EC2
if(codelb.ge.1.0) then
  codelb=0.
else if(codelb.le.-1.0) then
  codelb=pi
else
  DEL=ACOS(CODELB)
endif
AS=SC1*SC4
BS=SC1*SC3
AZI1=(AS-EC3)**2+(BS+EC4)**2+SC2*SC2-2.
AZI2=(AS-EC2*EC4)**2+(BS-EC2*EC3)**2+(SC2+EC1)**2-2.
IF(abs(AZI2)<1.0d-10) then
  AZIE=3.141592653-SIGN(1.570796326,AZI1)
else
  AZIE=ATAN2(AZI1,AZI2)
endif
del=del*rad
azis=azis*rad
azie=azie*rad
RETURN
END
