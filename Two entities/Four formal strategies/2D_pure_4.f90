! nm - number of cases evaluated
! ax - a vector of characteristics for entity 1 in each evaluated case
! bx - a vector of characteristics for entity 2 in each evaluated case
! pnash - a vector of the resulting NE probability values (calculated for pure strategies) in the individual cases evaluated
 
  subroutine d2pur4(nm,ax,bx,pnash)

	implicit real*8 (a-h,o-z)

	dimension avect(4,4),bvect(4,4)
	dimension inea(4),ineb(4)
	dimension ax(120),bx(120),axi(120),bxi(120),af(10),bf(10),pnash(120)

! matrix A dimension
     iar=4

! factor - eq.(2) * eq (9)
			nh=nm-1
			na=iar-1
			nd=na*2
			pfact=pbinom(nh,nd)*pbinom(nd,na)



 ! evaluated case 
	do 102 inasp=1,nm
	ipnas=0


! reindexing others cases
	in=1
	do 602 ir=1,inasp-1
	axi(in)=ax(ir)
	bxi(in)=bx(ir)
602	in=in+1
	do 702 ir=inasp+1,nm
	axi(in)=ax(ir)
	bxi(in)=bx(ir)
702	in=in+1
	

! selection (A-1)*N cases

	do 156 i=1,nm-1
	 do 156 j=i+1,nm-1
	  do 156 k=j+1,nm-1
	   do 156 l=k+1,nm-1
	    do 156 m=l+1,nm-1
	     do 156 n=m+1,nm-1

	   af(1)=axi(i)
	   af(2)=axi(j)
	   af(3)=axi(k)
	   af(4)=axi(l)
	   af(5)=axi(m)
	   af(6)=axi(n)


	   bf(1)=bxi(i)
	   bf(2)=bxi(j)
	   bf(3)=bxi(k)
	   bf(4)=bxi(l)
	   bf(5)=bxi(m)
	   bf(6)=bxi(n)



 ! combination into vectors - eq. (3)
	   avect(1,1)=af(1)
	   bvect(1,1)=bf(1)

	   do 502 ia2=2,6 

	   avect(1,2)=af(ia2)
	   bvect(1,2)=bf(ia2)

		do 502 ia3=ia2+1,6

	    avect(1,3)=af(ia3)
	    bvect(1,3)=bf(ia3)


	    do 402 ir=2,ia2-1
        avect(2,1)=af(ir)
	    bvect(2,1)=bf(ir)
		irr=ir
		goto 303
402	    continue
	    do 302 ir=ia2+1,ia3-1
        avect(2,1)=af(ir)
	    bvect(2,1)=bf(ir)
		irr=ir
		goto 303
302	    continue
	    do 202 ir=ia3+1,6
        avect(2,1)=af(ir)
	    bvect(2,1)=bf(ir)
		irr=ir
		goto 303
202	    continue
303 continue


 		ihf2=(abs(ia2+1-(irr+1))-(ia2+1-(irr+1)))/2
		ihf3=(abs(ia3+1-(irr+1))-(ia3+1-(irr+1)))/2


	    do 403 ir=irr+1,ia2-1
        avect(2,2)=af(ir)
	    bvect(2,2)=bf(ir)
		irr2=ir
		goto 304
403	    continue
	    do 405 ir=ia2+1+ihf2,ia3-1
        avect(2,2)=af(ir)
	    bvect(2,2)=bf(ir)
		irr2=ir
		goto 304
405	    continue
	    do 203 ir=ia3+1+ihf3,6
        avect(2,2)=af(ir)
	    bvect(2,2)=bf(ir)
		irr2=ir
		goto 304
203	    continue
304 continue


		iap=21-1-ia2-ia3-irr-irr2

	    avect(2,3)=af(iap)
	    bvect(2,3)=bf(iap)


 ! permutations of vectors - eq. (5) NE evaluation - eq (6),(7)
 		 do 305 iva=1,2

			do 306 ivv=1,3
			 inea(ivv)=0
			 if(ax(inasp).ge.avect(iva,ivv)) inea(ivv)=1
306			continue
			innea=inea(1)*inea(2)*inea(3)


		  do 308 ivb=1,iva-1

			do 307 ivv=1,3
			 ineb(ivv)=0
			 if(bx(inasp).ge.bvect(ivb,ivv)) ineb(ivv)=1
307			continue
			inneb=ineb(1)*ineb(2)*ineb(3)

308			ipnas=ipnas+innea*inneb			       


		   do 309 ivb=iva+1,2

			do 310 ivv=1,3
			 ineb(ivv)=0
			 if(bx(inasp).ge.bvect(ivb,ivv)) ineb(ivv)=1
310			continue
			inneb=ineb(1)*ineb(2)*ineb(3)

309			ipnas=ipnas+innea*inneb		   ! eq. (8),(10)

305      continue

502   continue 
156  continue

		pnash(inasp)=ipnas


102  continue

 	do 126 i=1,nm
 	pnash(i)=pnash(i)/pfact		     	  ! eq. (11)
126 continue 

     return
	 end
