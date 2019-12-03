! nm - number of cases evaluated
! ax - a vector of characteristics for entity 1 in each evaluated case
! bx - a vector of characteristics for entity 2 in each evaluated case
! cx - a vector of characteristics for entity 3 in each evaluated case
! pnash - a vector of the resulting NE probability values (calculated for pure strategies) in the individual cases evaluated

	subroutine  d3pur3(nm,ax,bx,cx,pnash)

	implicit real*8 (a-h,o-z)

	dimension avect(4,4),bvect(4,4),cvect(4,4),af1(10),bf1(10),cf1(10)
	dimension inea(4),ineb(4),inec(4),ivbd(4),ivcd(4)
	dimension ax(120),bx(120),cx(120),axi(120),bxi(120),cxi(120),af(10),bf(10),cf(10),pnash(120)


	nvd=1

! 3 entities
     inh=3

! matrix A dimension
     iar=3


! factor - eq.(2) * eq (9)
			nh=nm-1
			na=iar-1
			nd1=na*inh
			nd2=nd1-na
			pfact=pbinom(nh,nd1)*pbinom(nd1,na)*pbinom(nd2,na)

 
! evaluated case
 	do 102 inasp=1,nm
	ipnas=0

! reindexing others cases
	in=1
	do 602 ir=1,inasp-1
	axi(in)=ax(ir)
	bxi(in)=bx(ir)
	cxi(in)=cx(ir)
602	in=in+1
	do 702 ir=inasp+1,nm
	axi(in)=ax(ir)
	bxi(in)=bx(ir)
	cxi(in)=cx(ir)
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


	   cf(1)=cxi(i)
	   cf(2)=cxi(j)
	   cf(3)=cxi(k)
	   cf(4)=cxi(l)
	   cf(5)=cxi(m)
	   cf(6)=cxi(n)


	   avect(1,1)=af(1)
	   bvect(1,1)=bf(1)
	   cvect(1,1)=cf(1)

	   do 522 ia1=2,6 

	   avect(1,2)=af(ia1)
	   bvect(1,2)=bf(ia1)
	   cvect(1,2)=cf(ia1)

!  reindexation of remaining cases
		irn=0
	    do 422 ir=2,ia1-1
        irn=irn+1
		af1(irn)=af(ir)
	    bf1(irn)=bf(ir)
	    cf1(irn)=cf(ir)
422	    continue
	    do 322 ir=ia1+1,6
        irn=irn+1
		af1(irn)=af(ir)
	    bf1(irn)=bf(ir)
	    cf1(irn)=cf(ir)
322	    continue


 ! combination into vectors - eq. (3)
 	   avect(2,1)=af1(1)
	   bvect(2,1)=bf1(1)
   	   cvect(2,1)=cf1(1)

	   do 502 ia=2,4 

	   avect(2,2)=af1(ia)
	   bvect(2,2)=bf1(ia)
	   cvect(2,2)=cf1(ia)


	    do 402 ir=2,ia-1

        avect(3,1)=af1(ir)
	    bvect(3,1)=bf1(ir)
	    cvect(3,1)=cf1(ir)
		irr=ir
		goto 303
402	    continue
	    do 302 ir=ia+1,4
        avect(3,1)=af1(ir)
	    bvect(3,1)=bf1(ir)
	    cvect(3,1)=cf1(ir)
		irr=ir
		goto 303
302	    continue

303 continue

		iap=10-1-ia-irr

	    avect(3,2)=af1(iap)
	    bvect(3,2)=bf1(iap)
  	    cvect(3,2)=cf1(iap)

 
 ! permutations of vectors - eq. (5) NE evaluation - eq (6),(7)
 		 do 305 iva=1,inh

			do 306 ivv=1,2
			 inea(ivv)=0
			 if(ax(inasp).ge.avect(iva,ivv)) inea(ivv)=1
306			continue
			innea=inea(1)*inea(2)

		  id=1
		  do 408 iv=1,iva-1
			 ivbd(id)=iv
408			 id=id+1
		  do 409 iv=iva+1,3
			 ivbd(id)=iv
409			 id=id+1

		   do 305 idb=1,inh-1
			ivb=ivbd(idb)

		   	do 307 ivv=1,2
			 ineb(ivv)=0
			 if(bx(inasp).ge.bvect(ivb,ivv)) ineb(ivv)=1
307			continue
			inneb=ineb(1)*ineb(2)

		  id=1
		  do 406 iv=1,idb-1
			 ivcd(id)=iv
406			 id=id+1
		  do 407 iv=idb+1,3-1
			 ivcd(id)=iv
407			 id=id+1


		   do 305 idc=1,inh-2
		   
			ivc=ivcd(idc)

		   	do 309 ivv=1,2
			 inec(ivv)=0
			 if(cx(inasp).ge.cvect(ivc,ivv)) inec(ivv)=1
309			continue
			innec=inec(1)*inec(2)

			ipnas=ipnas+innea*inneb*innec	    ! eq. (8),(10)

305      continue

502   continue 
522  continue
156  continue

			pnash(inasp)=ipnas

102  continue

 	do 126 i=1,nm
 	pnash(i)=pnash(i)/pfact				    	! eq. (11)
126 continue 

     return
	 end

