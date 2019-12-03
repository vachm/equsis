! nm - number of cases evaluated
! ax - a vector of characteristics for entity 1 in each evaluated case
! bx - a vector of characteristics for entity 2 in each evaluated case
! cx - a vector of characteristics for entity 2 in each evaluated case
! slbxx - a vector of the resulting NE probability values (calculated for mixed strategies) for entity 1 in the individual cases evaluated
! slbyy - a vector of the resulting NE probability values (calculated for mixed strategies) for entity 2 in the individual cases evaluated
! slbzz - a vector of the resulting NE probability values (calculated for mixed strategies) for entity 3 in the individual cases evaluated

    subroutine d3mix2(nm,ax,bx,cx,slbxx,slbyy,slbzz)

  implicit real*8 (a-h,o-z)
	implicit integer*4 (i-n)

	dimension ax(120),bx(120),cx(120),axa(20),bxa(20),cxa(20),slbx(20),slby(20),slbz(20),slbxx(120),slbyy(120),slbzz(120)
	dimension nd(10),iid(10),jjd(10),jpd(10),kpd(10)

	zero=0.00000000000000

! 3 entities
     inh=3
	 ndim=inh

! matrix A dimension
     iar=2

! factor - eq.(12) * eq (13)	
			nhor=nm-1
			nmat=iar**ndim
			ndol=nmat-1
			pfact=pbinom(nhor,ndol)*pfc(nmat)/(pfc(iar)**ndim)

	nma=nmat

	nvd=1

! evaluated cases - 2 x 2 x 2 configuration 
	do 206 i1=1,nm
	do 206 i2=i1+1,nm 
	do 206 i3=i2+1,nm 
	do 206 i4=i3+1,nm 
	do 206 i5=i4+1,nm 
	do 206 i6=i5+1,nm 
	do 206 i7=i6+1,nm 
	do 206 i8=i7+1,nm 


! reindexing selected cases
	axa(1)=ax(i1)
 	axa(2)=ax(i2)
	axa(3)=ax(i3)
	axa(4)=ax(i4)
	axa(5)=ax(i5)
	axa(6)=ax(i6)
	axa(7)=ax(i7)
	axa(8)=ax(i8)

	bxa(1)=bx(i1)
 	bxa(2)=bx(i2)
	bxa(3)=bx(i3)
	bxa(4)=bx(i4)
	bxa(5)=bx(i5)
	bxa(6)=bx(i6)
	bxa(7)=bx(i7)
	bxa(8)=bx(i8)

	cxa(1)=cx(i1)
 	cxa(2)=cx(i2)
	cxa(3)=cx(i3)
	cxa(4)=cx(i4)
	cxa(5)=cx(i5)
	cxa(6)=cx(i6)
	cxa(7)=cx(i7)
	cxa(8)=cx(i8)


	do 562 it=1,nma
	 slbx(it)=zero
	 slby(it)=zero
	 slbz(it)=zero
562 continue
 

   i=1

! permutations - eq. (15), structure (14)
	do 102 l=2,nma

	in=1
	do 602 ir=2,l-1
	nd(in)=ir
602	in=in+1
	do 702 ir=l+1,nma
	nd(in)=ir
702	in=in+1
	
	 do	102 ni=1,nma-2
	 n=nd(ni)


	in=1
	do 603 ir=1,ni-1
	iid(in)=nd(ir)
603	in=in+1
	do 703 ir=ni+1,nma-2
	iid(in)=nd(ir)
703	in=in+1

	  do 102 iii=1,nma-3
	  ii=iid(iii)


	in=1
	do 604 ir=1,iii-1
	jjd(in)=iid(ir)
604	in=in+1
	do 704 ir=iii+1,nma-3
	jjd(in)=iid(ir)
704	in=in+1

	   do 102 jji=1,nma-4
	   jj=jjd(jji)


	in=1
	do 605 ir=1,jji-1
	jpd(in)=jjd(ir)
605	in=in+1
	do 705 ir=jji+1,nma-4
	jpd(in)=jjd(ir)
705	in=in+1

		  do 106 jp=1,3
			 j=jpd(jp)



	in=1
	do 606 ir=1,jp-1
	kpd(in)=jpd(ir)
606	in=in+1
	do 706 ir=jp+1,nma-5
	kpd(in)=jpd(ir)
706	in=in+1

		  do 106 kp=1,2
			 k=kpd(kp)

 			   m=36-i-j-k-l-n-ii-jj

  									  ihd=ihd+1


! NE evaluation - eq. (16)

  a111=axa(i)
  a121=axa(j)
  a211=axa(k)
  a221=axa(l)
  a112=axa(m)
  a122=axa(n)
  a212=axa(ii)
  a222=axa(jj)

  b111=bxa(i)
  b121=bxa(j)
  b211=bxa(k)
  b221=bxa(l)
  b112=bxa(m)
  b122=bxa(n)
  b212=bxa(ii)
  b222=bxa(jj)

  c111=cxa(i)
  c121=cxa(j)
  c211=cxa(k)
  c221=cxa(l)
  c112=cxa(m)
  c122=cxa(n)
  c212=cxa(ii)
  c222=cxa(jj)


  g1=c222-c212-c221+c211
  g2=-c222+c221
  g3=c222-c212-c122+c112-c221+c211+c121-c111
  g4=-c222+c122+c221-c121
  h1=a222-a122
  h2=a221-a121
  h3=a222-a212-a122+a112
  h4=a221-a211-a121+a111
  h5=b222-b212
  h6=b221-b211
  h7=b222-b212-b122+b112
  h8=b221-b211-b121+b111


  p11=-1
  p12=-1
  p21=-1
  p22=-1

  q11=-1
  q12=-1
  q21=-1
  q22=-1

  r11=-1
  r12=-1
  r21=-1
  r22=-1

  sqr=(g2**2*h3**2+2*g1*g2*h1*h3+g1**2*h1**2)*h8**2+((((-2*g2**2*h3)-2*g1*g2*h1)*h4-2*g1*g2*h2*h3-2*g1**2*h1*h2)*h7+((-2*g2*g4*h3**2)+((-2*g1*g4)-2*g2*g3)*h1*h3-2*g1*g3*h1**2)*h6+((2*g2*g4*h3+(4*g2*g3-2*g1*g4)*h1)*h4+(4*g1*g4-2*g2*g3)*h2*h3+2*g1*g3*h1*h2)*h5)*h8+(g2**2*h4**2+2*g1*g2*h2*h4+g1**2*h2**2)*h7**2+(((2*g2*g4*h3+(4*g1*g4-2*g2*g3)*h1)*h4+(4*g2*g3-2*g1*g4)*h2*h3+2*g1*g3*h1*h2)*h6+((-2*g2*g4*h4**2)+((-2*g1*g4)-2*g2*g3)*h2*h4-2*g1*g3*h2**2)*h5)*h7+(g4**2*h3**2+2*g3*g4*h1*h3+g3**2*h1**2)*h6**2+(((-2*g4**2*h3)-2*g3*g4*h1)*h4-2*g3*g4*h2*h3-2*g3**2*h1*h2)*h5*h6+(g4**2*h4**2+2*g3*g4*h2*h4+g3**2*h2**2)*h5**2
  dvr=(2*g2*h4-2*g2*h3+2*g1*h2-2*g1*h1)*h8+((-2*g2*h4)+2*g2*h3-2*g1*h2+2*g1*h1)*h7+((-2*g4*h4)+2*g4*h3-2*g3*h2+2*g3*h1)*h6+(2*g4*h4-2*g4*h3+2*g3*h2-2*g3*h1)*h5

  if(sqr.ge.zero.and.dvr.ne.zero) then
   
   r11=-(sqrt(sqr)+(g2*h3+g1*h1)*h8+(g2*h4-2*g2*h3+g1*h2-2*g1*h1)*h7+((-g4*h3)-g3*h1)*h6+((-g4*h4)+2*g4*h3-g3*h2+2*g3*h1)*h5)/dvr
   r12=(sqrt(sqr)+((-g2*h3)-g1*h1)*h8+((-g2*h4)+2*g2*h3-g1*h2+2*g1*h1)*h7+(g4*h3+g3*h1)*h6+(g4*h4-2*g4*h3+g3*h2-2*g3*h1)*h5)/dvr

  end if

  dvq11=h3*(1-r11)+h4*r11
  dvq12=h3*(1-r12)+h4*r12

  dvp11=h7*(1-r11)+h8*r11
  dvp12=h7*(1-r12)+h8*r12
 
  if(dvq11.ne.zero) q11=(h1*(1-r11)+h2*r11)/dvq11
  if(dvq12.ne.zero) q12=(h1*(1-r12)+h2*r12)/dvq12
 
  if(dvp11.ne.zero) p11=(h5*(1-r11)+h6*r11)/dvp11
  if(dvp12.ne.zero) p12=(h5*(1-r12)+h6*r12)/dvp12

  r21=1-r11
  r22=1-r12
  q21=1-q11
  q22=1-q12
  p21=1-p11
  p22=1-p12


  nash1=0
  nash2=0

  if(p11.ge.zero.and.q11.ge.zero.and.p21.ge.zero.and.q21.ge.zero.and.r11.ge.zero.and.r21.ge.zero) nash1=1
  if(p12.ge.zero.and.q12.ge.zero.and.p22.ge.zero.and.q22.ge.zero.and.r12.ge.zero.and.r22.ge.zero) nash2=1

   nashs=nash1+nash2

     if(nashs.eq.zero) nashs=1 

   nashf1=nash1/nashs
   nashf2=nash2/nashs


 		   slbx(i)=slbx(i)+(p11*nashf1+p12*nashf2)/4		     ! eq. (17)
 		   slbx(j)=slbx(j)+(p11*nashf1+p12*nashf2)/4
 		   slbx(m)=slbx(m)+(p11*nashf1+p12*nashf2)/4
 		   slbx(n)=slbx(n)+(p11*nashf1+p12*nashf2)/4

 		   slbx(k)=slbx(k)+(p21*nashf1+p22*nashf2)/4
 		   slbx(l)=slbx(l)+(p21*nashf1+p22*nashf2)/4
 		   slbx(ii)=slbx(ii)+(p21*nashf1+p22*nashf2)/4
 		   slbx(jj)=slbx(jj)+(p21*nashf1+p22*nashf2)/4


   		   slby(i)=slby(i)+(q11*nashf1+q12*nashf2)/4
   		   slby(k)=slby(k)+(q11*nashf1+q12*nashf2)/4
   		   slby(m)=slby(m)+(q11*nashf1+q12*nashf2)/4
   		   slby(ii)=slby(ii)+(q11*nashf1+q12*nashf2)/4

   		   slby(j)=slby(j)+(q21*nashf1+q22*nashf2)/4
   		   slby(l)=slby(l)+(q21*nashf1+q22*nashf2)/4
   		   slby(n)=slby(n)+(q21*nashf1+q22*nashf2)/4
   		   slby(jj)=slby(jj)+(q21*nashf1+q22*nashf2)/4


   		   slbz(i)=slbz(i)+(r11*nashf1+r12*nashf2)/4
   		   slbz(j)=slbz(j)+(r11*nashf1+r12*nashf2)/4
   		   slbz(k)=slbz(k)+(r11*nashf1+r12*nashf2)/4
   		   slbz(l)=slbz(l)+(r11*nashf1+r12*nashf2)/4

   		   slbz(m)=slbz(m)+(r21*nashf1+r22*nashf2)/4
   		   slbz(n)=slbz(n)+(r21*nashf1+r22*nashf2)/4
   		   slbz(ii)=slbz(ii)+(r21*nashf1+r22*nashf2)/4
   		   slbz(jj)=slbz(jj)+(r21*nashf1+r22*nashf2)/4

 
 106 continue
102  continue


	slbxx(i1)=slbxx(i1)+slbx(1)					 ! eq. (18)
	slbxx(i2)=slbxx(i2)+slbx(2)
	slbxx(i3)=slbxx(i3)+slbx(3)
	slbxx(i4)=slbxx(i4)+slbx(4)
	slbxx(i5)=slbxx(i5)+slbx(5)
	slbxx(i6)=slbxx(i6)+slbx(6)
	slbxx(i7)=slbxx(i7)+slbx(7)
	slbxx(i8)=slbxx(i8)+slbx(8)

	slbyy(i1)=slbyy(i1)+slby(1)
	slbyy(i2)=slbyy(i2)+slby(2)
	slbyy(i3)=slbyy(i3)+slby(3)
	slbyy(i4)=slbyy(i4)+slby(4)
	slbyy(i5)=slbyy(i5)+slby(5)
	slbyy(i6)=slbyy(i6)+slby(6)
	slbyy(i7)=slbyy(i7)+slby(7)
	slbyy(i8)=slbyy(i8)+slby(8)

	slbzz(i1)=slbzz(i1)+slbz(1)
	slbzz(i2)=slbzz(i2)+slbz(2)
	slbzz(i3)=slbzz(i3)+slbz(3)
	slbzz(i4)=slbzz(i4)+slbz(4)
	slbzz(i5)=slbzz(i5)+slbz(5)
	slbzz(i6)=slbzz(i6)+slbz(6)
	slbzz(i7)=slbzz(i7)+slbz(7)
	slbzz(i8)=slbzz(i8)+slbz(8)


206 continue

 	do 126 i=1,nm
 	  slbxx(i)=slbxx(i)/pfact					  ! eq. (21)
 	  slbyy(i)=slbyy(i)/pfact
 	  slbzz(i)=slbzz(i)/pfact
126 continue 

	return
	end
