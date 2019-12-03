! nm - number of cases evaluated
! ax - a vector of characteristics for entity 1 in each evaluated case
! bx - a vector of characteristics for entity 2 in each evaluated case
! slbxx - a vector of the resulting NE probability values (calculated for mixed strategies) for entity 1 in the individual cases evaluated
! slbyy - a vector of the resulting NE probability values (calculated for mixed strategies) for entity 2 in the individual cases evaluated


  subroutine d2mix3(nm,ax,bx,slbxx,slbyy)

  implicit real*8 (a-h,o-z)
	implicit integer*4 (i-n)

	dimension ax(120),bx(120),axa(20),bxa(20),slbx(20),slby(20),slbxx(120),slbyy(120)
	dimension jjd(20),kkd(20),ipd(20),iil(20)
	dimension nd(20)

	zero=0.000000000000

! matrix Amax dimension
	iar=3

! 2 entities
	ndim=2

! factor - eq.(12) * eq (13)
			nhor=nm-1
			nmat=iar**ndim
			ndol=nmat-1
			pfact=pbinom(nhor,ndol)*pfc(nmat)/(pfc(iar)**ndim)

 	nma=nmat


! evaluated cases - 3 x 3 matrix 
	do 206 i1=1,nm
	do 206 i2=i1+1,nm 
	do 206 i3=i2+1,nm 
	do 206 i4=i3+1,nm 
	do 206 i5=i4+1,nm 
	do 206 i6=i5+1,nm 
	do 206 i7=i6+1,nm 
	do 206 i8=i7+1,nm 
	do 206 i9=i8+1,nm 

 
! reindexing selected cases
	axa(1)=ax(i1)
 	axa(2)=ax(i2)
	axa(3)=ax(i3)
	axa(4)=ax(i4)
	axa(5)=ax(i5)
	axa(6)=ax(i6)
	axa(7)=ax(i7)
	axa(8)=ax(i8)
	axa(9)=ax(i9)

	bxa(1)=bx(i1)
 	bxa(2)=bx(i2)
	bxa(3)=bx(i3)
	bxa(4)=bx(i4)
	bxa(5)=bx(i5)
	bxa(6)=bx(i6)
	bxa(7)=bx(i7)
	bxa(8)=bx(i8)
	bxa(9)=bx(i9)


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
	
	 do	102 ni=1,nma-2
	 n=nd(ni)


	in=1
	do 603 ir=1,ni-1
	jjd(in)=nd(ir)
603	in=in+1
	do 703 ir=ni+1,nma-2
	jjd(in)=nd(ir)
703	in=in+1

	  do 102 jji=1,nma-3
	  jj=jjd(jji)


	in=1
	do 604 ir=1,jji-1
	kkd(in)=jjd(ir)
604	in=in+1
	do 704 ir=jji+1,nma-3
	kkd(in)=jjd(ir)
704	in=in+1

	   do 102 kki=1,nma-4
	   kk=kkd(kki)


	in=1
	do 605 ir=1,kki-1
	ipd(in)=kkd(ir)
605	in=in+1
	do 705 ir=kki+1,nma-4
	ipd(in)=kkd(ir)
705	in=in+1


		  do 102 jp=1,4
		   jjp=jp+1
		   do 102 kp=jjp,4

			 j=ipd(jp)
			 k=ipd(kp)

			iiil=1
			do 306 lp=1,jp-1
			 iil(iiil)=lp
			 iiil=iiil+1
306			continue 
			do 307 lp=jp+1,kp-1
			 iil(iiil)=lp
			 iiil=iiil+1
307			continue 
			do 308 lp=kp+1,4
			 iil(iiil)=lp
			 iiil=iiil+1
308			continue 


			 lind=iil(1)
			 iiind=iil(2)
			 
			 l=ipd(lind)
			 ii=ipd(iiind)
 
! NE evaluation - eq. (16)

		   p1=-1
		   p2=-1
		   p3=-1
		   q1=-1
		   q2=-1
		   q3=-1


		   cq1=(axa(m)-axa(j))*axa(kk)+(axa(k)-axa(n))*axa(jj)+axa(j)*axa(n)-axa(k)*axa(m)
		   cq2=-((axa(l)-axa(i))*axa(kk)+(axa(k)-axa(n))*axa(ii)+axa(i)*axa(n)-axa(k)*axa(l))

		   dq=(axa(m)-axa(l)-axa(j)+axa(i))*axa(kk)+(-axa(n)+axa(l)+axa(k)-axa(i))*axa(jj)+(axa(n)-axa(m)-axa(k)+axa(j))*axa(ii)+(axa(j)-axa(i))*axa(n)+(axa(i)-axa(k))*axa(m)+(axa(k)-axa(j))*axa(l)

		   if(dq.ne.zero) then
		   q1=cq1/dq
		   q2=cq2/dq
		   q3=1-q1-q2
				        end if
 
		   cp1=(bxa(m)-bxa(l))*bxa(kk)+(bxa(l)-bxa(n))*bxa(jj)+(bxa(n)-bxa(m))*bxa(ii)
		   cp2=-((bxa(j)-bxa(i))*bxa(kk)+(bxa(i)-bxa(k))*bxa(jj)+(bxa(k)-bxa(j))*bxa(ii))
		   
		   dp=(bxa(m)-bxa(l)-bxa(j)+bxa(i))*bxa(kk)+(-bxa(n)+bxa(l)+bxa(k)-bxa(i))*bxa(jj)+(bxa(n)-bxa(m)-bxa(k)+bxa(j))*bxa(ii)+(bxa(j)-bxa(i))*bxa(n)+(bxa(i)-bxa(k))*bxa(m)+(bxa(k)-bxa(j))*bxa(l)
			  
		   if(dp.ne.zero)	then
		   p1=cp1/dp
		   p2=cp2/dp
		   p3=1-p1-p2
						end if
	   		  
 		   if(p1.ge.zero.and.q1.ge.zero.and.p2.ge.zero.and.q2.ge.zero.and.p3.ge.zero.and.q3.ge.zero) then
 
		   slbx(i)=slbx(i)+p1/iar	 	! eq. (17)
		   slbx(j)=slbx(j)+p1/iar
		   slbx(k)=slbx(k)+p1/iar
		   slbx(l)=slbx(l)+p2/iar
		   slbx(m)=slbx(m)+p2/iar
		   slbx(n)=slbx(n)+p2/iar
		   slbx(ii)=slbx(ii)+p3/iar
		   slbx(jj)=slbx(jj)+p3/iar
		   slbx(kk)=slbx(kk)+p3/iar


		   slby(i)=slby(i)+q1/iar
		   slby(l)=slby(l)+q1/iar
		   slby(ii)=slby(ii)+q1/iar
		   slby(j)=slby(j)+q2/iar
		   slby(m)=slby(m)+q2/iar
		   slby(jj)=slby(jj)+q2/iar
		   slby(k)=slby(k)+q3/iar
		   slby(n)=slby(n)+q3/iar
		   slby(kk)=slby(kk)+q3/iar


                                       end if


102 continue

	slbxx(i1)=slbxx(i1)+slbx(1)				 ! eq. (18)
	slbxx(i2)=slbxx(i2)+slbx(2)
	slbxx(i3)=slbxx(i3)+slbx(3)
	slbxx(i4)=slbxx(i4)+slbx(4)
	slbxx(i5)=slbxx(i5)+slbx(5)
	slbxx(i6)=slbxx(i6)+slbx(6)
	slbxx(i7)=slbxx(i7)+slbx(7)
	slbxx(i8)=slbxx(i8)+slbx(8)
	slbxx(i9)=slbxx(i9)+slbx(9)

	slbyy(i1)=slbyy(i1)+slby(1)
	slbyy(i2)=slbyy(i2)+slby(2)
	slbyy(i3)=slbyy(i3)+slby(3)
	slbyy(i4)=slbyy(i4)+slby(4)
	slbyy(i5)=slbyy(i5)+slby(5)
	slbyy(i6)=slbyy(i6)+slby(6)
	slbyy(i7)=slbyy(i7)+slby(7)
	slbyy(i8)=slbyy(i8)+slby(8)
	slbyy(i9)=slbyy(i9)+slby(9)

206 continue 

 	do 126 i=1,nm
 	  slbxx(i)=slbxx(i)/pfact		      ! eq. (21)
 	  slbyy(i)=slbyy(i)/pfact
126 continue 

	return
	end


