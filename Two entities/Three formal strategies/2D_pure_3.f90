
! nm - number of cases evaluated
! ax - a vector of characteristics for entity 1 in each evaluated case
! bx - a vector of characteristics for entity 2 in each evaluated case
! pnash - a vector of the resulting NE probability values (calculated for pure strategies) in the individual cases evaluated

subroutine d2pur3(nm,ax,bx,pnash)

	implicit real*8 (a-h,o-z)

	dimension avect(4,4),bvect(4,4)
	dimension inea(4),ineb(4)
	dimension ax(120),bx(120),axi(120),bxi(120),af(10),bf(10),pnash(120)


! matrix A dimension
     iar=3

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

	   af(1)=axi(i)
	   af(2)=axi(j)
	   af(3)=axi(k)
	   af(4)=axi(l)


	   bf(1)=bxi(i)
	   bf(2)=bxi(j)
	   bf(3)=bxi(k)
	   bf(4)=bxi(l)



 ! permutations to vectors - eq. (3)
 	   avect(1,1)=af(1)
	   bvect(1,1)=bf(1)

	   do 502 ia=2,4 

	   avect(1,2)=af(ia)
	   bvect(1,2)=bf(ia)



	    do 402 ir=2,ia-1

        avect(2,1)=af(ir)
	    bvect(2,1)=bf(ir)
		irr=ir
		goto 303
402	    continue
	    do 302 ir=ia+1,4
        avect(2,1)=af(ir)
	    bvect(2,1)=bf(ir)
		irr=ir
		goto 303
302	    continue

303 continue

		iap=10-1-ia-irr

	    avect(2,2)=af(iap)
	    bvect(2,2)=bf(iap)

 
 ! permutations of vectors - eq. (5) NE evaluation - eq (6),(7)
		 do 305 iva=1,2

			do 306 ivv=1,2
			 inea(ivv)=0
			 if(ax(inasp).ge.avect(iva,ivv)) inea(ivv)=1
306			continue
			innea=inea(1)*inea(2)


		  do 308 ivb=1,iva-1

			do 307 ivv=1,2
			 ineb(ivv)=0
			 if(bx(inasp).ge.bvect(ivb,ivv)) ineb(ivv)=1
307			continue
			inneb=ineb(1)*ineb(2)

308			ipnas=ipnas+innea*inneb


		   do 309 ivb=iva+1,2

			do 310 ivv=1,2
			 ineb(ivv)=0
			 if(bx(inasp).ge.bvect(ivb,ivv)) ineb(ivv)=1
310			continue
			inneb=ineb(1)*ineb(2)

309			ipnas=ipnas+innea*inneb		 ! eq. (8),(10)

305      continue

502   continue 
156  continue

			pnash(inasp)=ipnas

102  continue

 	do 126 i=1,nm
 	pnash(i)=pnash(i)/pfact		  ! eq. (11)
126 continue 

	return
	end

