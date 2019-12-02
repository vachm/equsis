! nm - number of cases evaluated
! ax - a vector of characteristics for entity 1 in each evaluated case
! bx - a vector of characteristics for entity 2 in each evaluated case
! pnash - a vector of the resulting NE probability values (calculated for pure strategies) in the individual cases evaluated


subroutine d2pur2(nm,ax,bx,pnash)

	implicit real*8 (a-h,o-z)

	dimension avect(4,4),bvect(4,4)
	dimension inea(4),ineb(4)
	dimension ax(120),bx(120),axi(120),bxi(120),af(5),bf(5),pnash(120)


! matrix A dimension
    iar=2

! factor - eq.(2) * eq (9)
 			nh=nm-1
			na=iar-1
			nd=na*iar
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

	   af(1)=axi(i)
	   af(2)=axi(j)

	   bf(1)=bxi(i)
	   bf(2)=bxi(j)


 ! permutations to vectors - eq. (3)
	   avect(1,1)=af(1)
	   bvect(1,1)=bf(1)

       avect(2,1)=af(2)
       bvect(2,1)=bf(2)


 
 ! permutations of vectors - eq. (5) NE evaluation - eq (6),(7)
		 do 305 iva=1,2

			 inea(1)=0
			 if(ax(inasp).ge.avect(iva,1)) inea(1)=1
			 innea=inea(1)


		  do 308 ivb=1,iva-1

			 ineb(1)=0
			 if(bx(inasp).ge.bvect(ivb,1)) ineb(1)=1
			inneb=ineb(1)

308			ipnas=ipnas+innea*inneb


		   do 309 ivb=iva+1,2

			 ineb(1)=0
			 if(bx(inasp).ge.bvect(ivb,1)) ineb(1)=1
			inneb=ineb(1)

309			ipnas=ipnas+innea*inneb	   ! eq. (8),(10)

305      continue

156  continue

			pnash(inasp)=ipnas

102  continue

 	do 126 i=1,nm
 	pnash(i)=pnash(i)/pfact		! eq. (11)
126 continue 

	return
	end

