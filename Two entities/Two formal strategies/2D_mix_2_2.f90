
! nm - number of cases evaluated
! ax - a vector of characteristics for entity (a) in each evaluated case
! bx - a vector of characteristics for entity (b) in each evaluated case
! slbxx - a vector of the resulting NE probability values (calculated for mixed strategies) for entity 1 in the individual cases evaluated
! slbyy - a vector of the resulting NE probability values (calculated for mixed strategies) for entity 2 in the individual cases evaluated


  subroutine d2mix2(nm,ax,bx,slbxx,slbyy)

 	implicit real*8 (a-h,o-z)
	implicit integer*4 (i-n)

	dimension ax(120),bx(120),axa(10),bxa(10),slbx(10),slby(10),slbxx(120),slbyy(120)
	dimension nd(10)


	zero=0.000000000000

! matrix A dimension
	iar=2
! 2 entities
	ndim=2

! factor - eq.(12) * eq (13)
			nhor=nm-1
			nmat=iar**ndim
			ndol=nmat-1
			pfact=pbinom(nhor,ndol)*pfc(nmat)/(pfc(iar)**ndim)



	nma=nmat

! evaluated cases - 2 x 2 matrix 
	do 206 i1=1,nm 
	do 206 i2=i1+1,nm 
	do 206 i3=i2+1,nm 
	do 206 i4=i3+1,nm 

! reindexing selected cases
	axa(1)=ax(i1)
 	axa(2)=ax(i2)
	axa(3)=ax(i3)
	axa(4)=ax(i4)

	bxa(1)=bx(i1)
 	bxa(2)=bx(i2)
	bxa(3)=bx(i3)
	bxa(4)=bx(i4)


	do 562 it=1,nma
	 slbx(it)=zero
	 slby(it)=zero
562 continue
 
   i=1
  

! permutations - eq. (15), structure (14)
	do 102 m=2,nma
	 
 	in=1
	do 602 ir=2,m-1
	nd(in)=ir
602	in=in+1
	do 702 ir=m+1,nma
	nd(in)=ir
702	in=in+1
	


		  do 106 jp=1,2

			 j=nd(jp)

			 l=9-m-j

					   ihd=ihd+1

 ! NE evaluation - eq. (16)

		   q1=eq1(axa(i),axa(j),axa(l),axa(m))
		   q2=1-q1

		   p1=ep1(bxa(i),bxa(j),bxa(l),bxa(m))
		   p2=1-p1
	   		  
 		   if(p1.ge.zero.and.q1.ge.zero.and.p2.ge.zero.and.q2.ge.zero) then
 
		   slbx(i)=slbx(i)+p1/ndim	! eq. (17)
		   slbx(j)=slbx(j)+p1/ndim
		   slbx(l)=slbx(l)+p2/ndim
		   slbx(m)=slbx(m)+p2/ndim

		   slby(i)=slby(i)+q1/ndim
		   slby(l)=slby(l)+q1/ndim
		   slby(j)=slby(j)+q2/ndim
		   slby(m)=slby(m)+q2/ndim

                                       end if

106 continue
102 continue



	slbxx(i1)=slbxx(i1)+slbx(1)		 ! eq. (18)
	slbxx(i2)=slbxx(i2)+slbx(2)
	slbxx(i3)=slbxx(i3)+slbx(3)
	slbxx(i4)=slbxx(i4)+slbx(4)

	slbyy(i1)=slbyy(i1)+slby(1)
	slbyy(i2)=slbyy(i2)+slby(2)
	slbyy(i3)=slbyy(i3)+slby(3)
	slbyy(i4)=slbyy(i4)+slby(4)


206 continue 

 	do 126 i=1,nm
 	  slbxx(i)=slbxx(i)/pfact	    ! eq. (21)
 	  slbyy(i)=slbyy(i)/pfact
126 continue 

   return
	  end





function eq1(t11,t12,t21,t22)
 	implicit real*8 (a-h,o-z)

	zero=0.00000000000000
    q1=-1
	
		   cq1=t22-t12
		   dq=t11-t12-t21+t22

		   if(dq.ne.zero) q1=cq1/dq
		   eq1=q1
	return
	end	
				 

function ep1(v11,v12,v21,v22)
 	implicit real*8 (a-h,o-z)

	zero=0.00000000000000
    p1=-1

		   cp1=v22-v21
		   dp=v22-v21-v12+v11
			  
		   if(dp.ne.zero) p1=cp1/dp
		   ep1=p1
 	return
	end	


