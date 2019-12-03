
! nm - number of cases evaluated
! ax - a vector of characteristics for entity 1 in each evaluated case
! bx - a vector of characteristics for entity 2 in each evaluated case
! cx - a vector of characteristics for entity 3 in each evaluated case
! dx - a vector of characteristics for entity 4 in each evaluated case
! pnash - a vector of the resulting NE probability values (calculated for pure strategies) in the individual cases evaluated

	subroutine  d3pur2(nm,ax,bx,cx,dx,pnash)
	implicit real*8 (a-h,o-z)

 	dimension avect(6,6),bvect(6,6),cvect(6,6),dvect(6,6)
	dimension inea(6),ineb(6),inec(6),ined(6),ivbd(6),ivcd(6),ivdd(6),af(10),bf(10),cf(10),df(10)
  dimension ax(120),bx(120),cx(120),dx(120),axi(120),bxi(120),cxi(120),dxi(120),pnash(120)

  nvd=1
  
! 4 entities
     inh=4
     
! matrix A dimension
     iar=2

! factor - eq.(2) * eq (9)
      nh=nm-1
			na=iar-1
			nd1=na*inh
			nd2=nd1-na
			nd3=nd2-na
			pfact=pbinom(nh,nd1)*pbinom(nd1,na)*pbinom(nd2,na)*pbinom(nd3,na)


! evaluated case 
	do 212 inasp=1,nm
	ipnas=0


! reindexing others cases
	in=1
	do 802 ir=1,inasp-1
	axi(in)=ax(ir)
	bxi(in)=bx(ir)
	cxi(in)=cx(ir)
	dxi(in)=dx(ir)
802	in=in+1

	do 902 ir=inasp+1,nm
	axi(in)=ax(ir)
	bxi(in)=bx(ir)
	cxi(in)=cx(ir)
	dxi(in)=dx(ir)
902	in=in+1
	

! selection (A-1)*N cases

	do 256 i=1,nm-1
	 do 256 j=i+1,nm-1
	  do 256 k=j+1,nm-1
	   do 256 l=k+1,nm-1

	   af(1)=axi(i)
	   af(2)=axi(j)
	   af(3)=axi(k)
	   af(4)=axi(l)

	   bf(1)=bxi(i)
	   bf(2)=bxi(j)
	   bf(3)=bxi(k)
	   bf(4)=bxi(l)

	   cf(1)=cxi(i)
	   cf(2)=cxi(j)
	   cf(3)=cxi(k)
	   cf(4)=cxi(l)

	   df(1)=dxi(i)
	   df(2)=dxi(j)
	   df(3)=dxi(k)
	   df(4)=dxi(l)

! combination into vectors - eq. (3)
     avect(1,1)=af(1)
	   bvect(1,1)=bf(1)
	   cvect(1,1)=cf(1)
	   dvect(1,1)=df(1)

	   avect(2,1)=af(2)
	   bvect(2,1)=bf(2)
	   cvect(2,1)=cf(2)
	   dvect(2,1)=df(2)

	   avect(3,1)=af(3)
	   bvect(3,1)=bf(3)
	   cvect(3,1)=cf(3)
	   dvect(3,1)=df(3)

	   avect(4,1)=af(4)
	   bvect(4,1)=bf(4)
	   cvect(4,1)=cf(4)
	   dvect(4,1)=df(4)

 ! permutations of vectors - eq. (5) NE evaluation - eq (6),(7)
		 do 415 iva=1,inh

			do 416 ivv=1,1
			 inea(ivv)=0
			 if(ax(inasp).ge.avect(iva,ivv)) inea(ivv)=1
416			continue
			innea=inea(1)


		  id=1
		  do 408 iv=1,iva-1
			 ivbd(id)=iv
408			 id=id+1
		  do 409 iv=iva+1,inh
			 ivbd(id)=iv
409			 id=id+1


		   do 415 idb=1,inh-1
		   
			ivb=ivbd(idb)

		   	do 327 ivv=1,1
			 ineb(ivv)=0
			 if(bx(inasp).ge.bvect(ivb,ivv)) ineb(ivv)=1
327			continue
			inneb=ineb(1)


		  id=1
		  do 406 iv=1,idb-1
			 ivcd(id)=iv
406			 id=id+1
		  do 407 iv=idb+1,inh-1
			 ivcd(id)=iv
407			 id=id+1



		   do 415 idc=1,inh-2
		   
			ivc=ivcd(idc)

		   	do 324 ivv=1,1
			 inec(ivv)=0
			 if(cx(inasp).ge.cvect(ivc,ivv)) inec(ivv)=1
324			continue
			innec=inec(1)


		  id=1
		  do 426 iv=1,idc-1
			 ivdd(id)=iv
426			 id=id+1
		  do 427 iv=idc+1,inh-2
			 ivdd(id)=iv
427			 id=id+1



		   do 415 idd=1,inh-3
		   
			ivd=ivdd(idd)

		   	do 329 ivv=1,1
			 ined(ivv)=0
			 if(dx(inasp).ge.dvect(ivd,ivv)) ined(ivv)=1
329			continue
			inned=ined(1)


			ipnas=ipnas+innea*inneb*innec*inned              ! eq. (8),(10)

415      continue


256  continue
			pnash(inasp)=ipnas


212  continue




 	do 646 i=1,nm
 	pnash(i)=pnash(i)/pfact             	! eq. (11)
646 continue 

	 return
	 end

