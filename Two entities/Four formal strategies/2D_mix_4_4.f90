    subroutine d2mix4(nm,ax,bx,slbxx,slbyy)

   	implicit real*8 (a-h,o-z)
	implicit integer*4 (i-n)

	dimension ax(120),bx(120),axa(20),bxa(20),slbx(20),slby(20),slbxx(120),slbyy(120)
	dimension jjd(20),ipd(20),iil(20),md(20),ld(20),iid(20),jd(20),id(20),kd(20),nd(20)

	zero=0.00000000000000

! matrix Amax dimension
	iar=4

! 2 entities
	ndim=2

! factor - eq.(12) * eq (13)
			nhor=nm-1
			nmat=iar**ndim
			ndol=nmat-1
			pfact=pbinom(nhor,ndol)*pfc(nmat)/(pfc(iar)**ndim)

	nma=nmat


! evaluated cases - 4 x 4 matrix 

 	do 206 i1=1,nm 
	do 206 i2=i1+1,nm 
	do 206 i3=i2+1,nm 
	do 206 i4=i3+1,nm 
	do 206 i5=i4+1,nm 
	do 206 i6=i5+1,nm 
	do 206 i7=i6+1,nm 
	do 206 i8=i7+1,nm 
	do 206 i9=i8+1,nm 
	do 206 i10=i9+1,nm 
	do 206 i11=i10+1,nm 
	do 206 i12=i11+1,nm 
	do 206 i13=i12+1,nm 
	do 206 i14=i13+1,nm 
	do 206 i15=i14+1,nm 
	do 206 i16=i15+1,nm 

 
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
	axa(10)=ax(i10)
	axa(11)=ax(i11)
	axa(12)=ax(i12)
	axa(13)=ax(i13)
	axa(14)=ax(i14)
	axa(15)=ax(i15)
	axa(16)=ax(i16)

	bxa(1)=bx(i1)
 	bxa(2)=bx(i2)
	bxa(3)=bx(i3)
	bxa(4)=bx(i4)
	bxa(5)=bx(i5)
	bxa(6)=bx(i6)
	bxa(7)=bx(i7)
	bxa(8)=bx(i8)
	bxa(9)=bx(i9)
	bxa(10)=bx(i10)
	bxa(11)=bx(i11)
	bxa(12)=bx(i12)
	bxa(13)=bx(i13)
	bxa(14)=bx(i14)
	bxa(15)=bx(i15)
	bxa(16)=bx(i16)



	do 562 it=1,nma
	 slbx(it)=zero
	 slby(it)=zero
562 continue


 
   ix0=1

! permutations - eq. (15), structure (14)
	do 102 i=2,nma
!		             write(*,*) i
	in=1
	do 602 ir=2,i-1
	id(in)=ir
602	in=in+1
	do 702 ir=i+1,nma
	id(in)=ir
702	in=in+1
	
	 do	102 ini=1,nma-2
	 j=id(ini)


	in=1
	do 603 ir=1,ini-1
	jd(in)=id(ir)
603	in=in+1
	do 703 ir=ini+1,nma-2
	jd(in)=id(ir)
703	in=in+1

	  do 102 ji=1,nma-3
	  k=jd(ji)

	in=1
	do 604 ir=1,ji-1
	kd(in)=jd(ir)
604	in=in+1
	do 704 ir=ji+1,nma-3
	kd(in)=jd(ir)
704	in=in+1

	   do 102 ki=1,nma-4
	   l=kd(ki)


	in=1
	do 605 ir=1,ki-1
	ld(in)=kd(ir)
605	in=in+1
	do 705 ir=ki+1,nma-4
	ld(in)=kd(ir)
705	in=in+1

	   do 102 li=1,nma-5
	   m=ld(li)


	in=1
	do 606 ir=1,li-1
	md(in)=ld(ir)
606	in=in+1
	do 706 ir=li+1,nma-5
	md(in)=ld(ir)
706	in=in+1

	   do 102 mi=1,nma-6
	   n=md(mi)


	in=1
	do 607 ir=1,mi-1
	nd(in)=md(ir)
607	in=in+1
	do 707 ir=mi+1,nma-6
	nd(in)=md(ir)
707	in=in+1

	   do 102 ni=1,nma-7
	   ii=nd(ni)


	in=1
	do 608 ir=1,ni-1
	iid(in)=nd(ir)
608	in=in+1
	do 708 ir=ni+1,nma-7
	iid(in)=nd(ir)
708	in=in+1

	   do 102 iii=1,nma-8
	   jj=iid(iii)


	in=1
	do 609 ir=1,iii-1
	jjd(in)=iid(ir)
609	in=in+1
	do 709 ir=iii+1,nma-8
	jjd(in)=iid(ir)
709	in=in+1

	   do 102 jji=1,nma-9
	   kk=jjd(jji)

										 	  ihd=ihd+1

									!		  write(*,*) ihd

	in=1
	do 610 ir=1,jji-1
	ipd(in)=jjd(ir)
610	in=in+1
	do 710 ir=jji+1,nma-9
	ipd(in)=jjd(ir)
710	in=in+1


		  do 102 ix1p=1,6
		   do 102 ix2p=ix1p+1,6
		    do 102 ix3p=ix2p+1,6

			 ix1=ipd(ix1p)
			 ix2=ipd(ix2p)
			 ix3=ipd(ix3p)

			iiil=1
			do 306 lp=1,ix1p-1
			 iil(iiil)=lp
			 iiil=iiil+1
306			continue 
			do 307 lp=ix1p+1,ix2p-1
			 iil(iiil)=lp
			 iiil=iiil+1
307			continue 
			do 308 lp=ix2p+1,ix3p-1
			 iil(iiil)=lp
			 iiil=iiil+1
308			continue 
			do 309 lp=ix3p+1,6
			 iil(iiil)=lp
			 iiil=iiil+1
309			continue 


			 il4=iil(1)
			 il5=iil(2)
			 il6=iil(3)
			 
			 ix4=ipd(il4)
			 ix5=ipd(il5)
			 ix6=ipd(il6)


										  xihd=xihd+1

! NE evaluation - eq. (16)

		   a11=axa(i)
		   a12=axa(j)
		   a13=axa(k)
		   a14=axa(ix1)

		   a21=axa(l)
		   a22=axa(m)
		   a23=axa(n)
		   a24=axa(ix2)

		   a31=axa(ii)
		   a32=axa(kk)
		   a33=axa(jj)
		   a34=axa(ix3)
		   
		   a41=axa(ix4)
		   a42=axa(ix5)
		   a43=axa(ix6)
		   a44=axa(ix0)
		   
		   b11=bxa(i)
		   b12=bxa(j)
		   b13=bxa(k)
		   b14=bxa(ix1)

		   b21=bxa(l)
		   b22=bxa(m)
		   b23=bxa(n)
		   b24=bxa(ix2)

		   b31=bxa(ii)
		   b32=bxa(kk)
		   b33=bxa(jj)
		   b34=bxa(ix3)
		   
		   b41=bxa(ix4)
		   b42=bxa(ix5)
		   b43=bxa(ix6)
		   b44=bxa(ix0)
		   
		   p1=-1
		   p2=-1
		   p3=-1
		   p4=-1
		   q1=-1
		   q2=-1
		   q3=-1
		   q4=-1



	   dq=(a21*(a12*(a44-a43-a34+a33)+a33*(a44-a42)+a32*(a43-a44)+a13*((-a44)+a42+a34-a32)+a14*(a43-a42-a33+a32)+a34*(a42-a43))+a11*(a32*(a44-a43)+a33*(a42-a44)+a34*(a43-a42))+a22*(a31*(a44-a43)+a13*(a44-a41-a34+a31)+a11*((-a44)+a43+a34-a33)+a33*(a41-a44)+a34*(a43-a41)+a14*((-a43)+a41+a33-a31))+a23*(a11*(a44-a42-a34+a32)+a32*(a44-a41)+a31*(a42-a44)+a12*((-a44)+a41+a34-a31)+a14*(a42-a41-a32+a31)+a34*(a41-a42))+a13*(a31*(a44-a42)+a32*(a41-a44)+a34*(a42-a41))+a12*(a33*(a44-a41)+a31*(a43-a44)+a34*(a41-a43))+a24*(a31*(a43-a42)+a12*(a43-a41-a33+a31)+a11*((-a43)+a42+a33-a32)+a32*(a41-a43)+a33*(a42-a41)+a13*((-a42)+a41+a32-a31))+a14*(a32*(a43-a41)+a31*(a42-a43)+a33*(a41-a42)))
	   cq1=(a22*(a13*(a44-a34)-a33*a44+a34*a43+a14*(a33-a43))+a12*(a33*a44-a34*a43)+a23*(a32*a44+a12*(a34-a44)+a14*(a42-a32)-a34*a42)+a13*(a34*a42-a32*a44)+a24*(a12*(a43-a33)-a32*a43+a33*a42+a13*(a32-a42))+a14*(a32*a43-a33*a42))
	   cq2=-(a21*(a13*(a44-a34)-a33*a44+a34*a43+a14*(a33-a43))+a11*(a33*a44-a34*a43)+a23*(a31*a44+a11*(a34-a44)+a14*(a41-a31)-a34*a41)+a13*(a34*a41-a31*a44)+a24*(a11*(a43-a33)-a31*a43+a33*a41+a13*(a31-a41))+a14*(a31*a43-a33*a41))
	   cq3=(a21*(a12*(a44-a34)-a32*a44+a34*a42+a14*(a32-a42))+a11*(a32*a44-a34*a42)+a22*(a31*a44+a11*(a34-a44)+a14*(a41-a31)-a34*a41)+a12*(a34*a41-a31*a44)+a24*(a11*(a42-a32)-a31*a42+a32*a41+a12*(a31-a41))+a14*(a31*a42-a32*a41))


 	   dp=(b21*(b12*(b44-b43-b34+b33)+b33*(b44-b14)+b13*(b34-b44)+(b14-b34)*b43)+b32*(b11*(b44-b43-b24+b23)+b23*(b44-b14)+b21*((-b44)+b43+b14-b13)+b13*(b24-b44)+(b14-b24)*b43)+b31*(b22*(b44-b43-b14+b13)+b13*(b44-b24)+b12*((-b44)+b43+b24-b23)+b23*(b14-b44)+(b24-b14)*b43)+b11*(b23*(b44-b34)+b33*(b24-b44)+(b34-b24)*b43)+b22*(b13*(b44-b34)+b11*((-b44)+b43+b34-b33)+b33*(b14-b44)+(b34-b14)*b43)+b12*(b33*(b44-b24)+b23*(b34-b44)+(b24-b34)*b43)+(b21*(b34-b33-b14+b13)+b13*(b34-b24)+b11*((-b34)+b33+b24-b23)+b23*(b14-b34)+(b24-b14)*b33+((-b24)+b23+b14-b13)*b31)*b42+(b12*(b34-b33-b24+b23)+b23*(b34-b14)+b22*((-b34)+b33+b14-b13)+b13*(b24-b34)+(b14-b24)*b33+(b24-b23-b14+b13)*b32)*b41)
	   cp1=(b31*(b22*(b44-b43)-b23*b44+b24*b43)+b21*(b33*b44-b34*b43)+b22*(b34*b43-b33*b44)+b32*(b23*b44+b21*(b43-b44)-b24*b43)+(b21*(b34-b33)-b23*b34+b24*b33+(b23-b24)*b31)*b42+(b23*b34+b22*(b33-b34)-b24*b33+(b24-b23)*b32)*b41)
	   cp2=-(b31*(b12*(b44-b43)-b13*b44+b14*b43)+b11*(b33*b44-b34*b43)+b12*(b34*b43-b33*b44)+b32*(b13*b44+b11*(b43-b44)-b14*b43)+(b11*(b34-b33)-b13*b34+b14*b33+(b13-b14)*b31)*b42+(b13*b34+b12*(b33-b34)-b14*b33+(b14-b13)*b32)*b41)
	   cp3=(b21*(b12*(b44-b43)-b13*b44+b14*b43)+b11*(b23*b44-b24*b43)+b12*(b24*b43-b23*b44)+b22*(b13*b44+b11*(b43-b44)-b14*b43)+(b11*(b24-b23)-b13*b24+b14*b23+(b13-b14)*b21)*b42+(b13*b24+b12*(b23-b24)-b14*b23+(b14-b13)*b22)*b41)

	   
	       if(dq.ne.zero) then

		   q1=cq1/dq
		   q2=cq2/dq
		   q3=cq3/dq
		   q4=1-q1-q2-q3
				        end if
			  
		   if(dp.ne.zero) then

		   p1=cp1/dp
		   p2=cp2/dp
		   p3=cp3/dp
		   p4=1-p1-p2-p3
						end if

	   		  
 		   if(p1.ge.zero.and.q1.ge.zero.and.p2.ge.zero.and.q2.ge.zero.and.p3.ge.zero.and.q3.ge.zero.and.p4.ge.zero.and.q4.ge.zero) then
 
		   slbx(i)=slbx(i)+p1/iar					  ! eq. (17)
		   slbx(j)=slbx(j)+p1/iar
		   slbx(k)=slbx(k)+p1/iar
   		   slbx(ix1)=slbx(ix1)+p1/iar

		   slbx(l)=slbx(l)+p2/iar
		   slbx(m)=slbx(m)+p2/iar
		   slbx(n)=slbx(n)+p2/iar
   		   slbx(ix2)=slbx(ix2)+p2/iar

		   slbx(ii)=slbx(ii)+p3/iar
		   slbx(jj)=slbx(jj)+p3/iar
		   slbx(kk)=slbx(kk)+p3/iar
   		   slbx(ix3)=slbx(ix3)+p3/iar

		   slbx(ix4)=slbx(ix4)+p4/iar
		   slbx(ix5)=slbx(ix5)+p4/iar
		   slbx(ix6)=slbx(ix6)+p4/iar
		   slbx(ix0)=slbx(ix0)+p4/iar


		   slby(i)=slby(i)+q1/iar
		   slby(l)=slby(l)+q1/iar
		   slby(ii)=slby(ii)+q1/iar
		   slby(ix4)=slby(ix4)+q1/iar

		   slby(j)=slby(j)+q2/iar
		   slby(m)=slby(m)+q2/iar
		   slby(jj)=slby(jj)+q2/iar
		   slby(ix5)=slby(ix5)+q2/iar

		   slby(k)=slby(k)+q3/iar
		   slby(n)=slby(n)+q3/iar
		   slby(kk)=slby(kk)+q3/iar
		   slby(ix6)=slby(ix6)+q3/iar

		   slby(ix1)=slby(ix1)+q4/iar
		   slby(ix2)=slby(ix2)+q4/iar
		   slby(ix3)=slby(ix3)+q4/iar
		   slby(ix0)=slby(ix0)+q4/iar

                                       end if




102 continue

 
	slbxx(i1)=slbxx(i1)+slbx(1)					  	 ! eq. (18)
	slbxx(i2)=slbxx(i2)+slbx(2)
	slbxx(i3)=slbxx(i3)+slbx(3)
	slbxx(i4)=slbxx(i4)+slbx(4)
	slbxx(i5)=slbxx(i5)+slbx(5)
	slbxx(i6)=slbxx(i6)+slbx(6)
	slbxx(i7)=slbxx(i7)+slbx(7)
	slbxx(i8)=slbxx(i8)+slbx(8)
	slbxx(i9)=slbxx(i9)+slbx(9)
	slbxx(i10)=slbxx(i10)+slbx(10)
	slbxx(i11)=slbxx(i11)+slbx(11)
	slbxx(i12)=slbxx(i12)+slbx(12)
	slbxx(i13)=slbxx(i13)+slbx(13)
	slbxx(i14)=slbxx(i14)+slbx(14)
	slbxx(i15)=slbxx(i15)+slbx(15)
	slbxx(i16)=slbxx(i16)+slbx(16)

	slbyy(i1)=slbyy(i1)+slby(1)
	slbyy(i2)=slbyy(i2)+slby(2)
	slbyy(i3)=slbyy(i3)+slby(3)
	slbyy(i4)=slbyy(i4)+slby(4)
	slbyy(i5)=slbyy(i5)+slby(5)
	slbyy(i6)=slbyy(i6)+slby(6)
	slbyy(i7)=slbyy(i7)+slby(7)
	slbyy(i8)=slbyy(i8)+slby(8)
	slbyy(i9)=slbyy(i9)+slby(9)
	slbyy(i10)=slbyy(i10)+slby(10)
	slbyy(i11)=slbyy(i11)+slby(11)
	slbyy(i12)=slbyy(i12)+slby(12)
	slbyy(i13)=slbyy(i13)+slby(13)
	slbyy(i14)=slbyy(i14)+slby(14)
	slbyy(i15)=slbyy(i15)+slby(15)
	slbyy(i16)=slbyy(i16)+slby(16)


206 continue

 	do 126 i=1,nm
 	  slbxx(i)=slbxx(i)/pfact				   ! eq. (21)
 	  slbyy(i)=slbyy(i)/pfact
126 continue 

    return
    end
