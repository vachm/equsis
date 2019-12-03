! nm - number of cases evaluated
! ax - a vector of characteristics for entity 1 in each evaluated case
! bx - a vector of characteristics for entity 2 in each evaluated case
! slbxx - a vector of the resulting NE probability values (calculated for mixed strategies) for entity 1 in the individual cases evaluated
! slbyy - a vector of the resulting NE probability values (calculated for mixed strategies) for entity 2 in the individual cases evaluated


  subroutine d2mix32(nm,ax,bx,slbxx,slbyy)

 	implicit real*8 (a-h,o-z)
	implicit integer*4 (i-n)

	dimension ax(120),bx(120),slbx2(120),slby2(120),slbxx(120),slbyy(120)
	dimension nd(10)
  
  call d2mix2(nm,ax,bx,slbx2,slby2)   ! 2 x 2 mix calculating
  
! matrix Amax dimension
	nsq=4
! matrix A dimension
	nmt=3
	
! factor - eq.(19) * eq (20)
   nsqs=nsq**2
	 nmts=nmt**2
	 fact3_4=(pbinom(nsq,nmt)**2)*(pbinom(nsq,nmt)**2)/pbinom(nsqs,nmts)
	 	
 	 do 989 i=1,nm
	     slbxx(i)=slbx2(i)*fact3_4
989	   slbyy(i)=slby2(i)*fact3_4

  return
  end
