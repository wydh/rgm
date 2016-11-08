*deck main
	program main
      implicit double precision (a-h,o-z)
      parameter(imt=62,jmt=142,imm=imt-1,imm1=imm-1,jmm=jmt-1,
     *     jmm1=jmm-1,kmt=120)
	common /vel/ u(imt,jmt),v(imt,jmt),ua(imt,jmt),va(imt,jmt),
     *  ub(imt,jmt),vb(imt,jmt),h(imt,jmt),ha(imt,jmt),hb(imt,jmt),
     *  hua(imt,jmt),hu(imt,jmt),hub(imt,jmt),f(jmt),gp(imt,jmt),
     *  gp0(imt,jmt),wsx(imt,jmt),tau(jmt),taud(jmt),dh0(imt,jmt)
	common /sca/acor,dxt,dxtr,dtts,dtuv,c2dtts,c2dtuv,drag,gpri,
     *         mxp0,mxp,mixeb,nfirst,nlast,iwrite,dzi,p_yr,year1,year2
	dimension psi(imt,jmt),eta(imt,jmt),d(kmt),hd(kmt)
     *            ,hct(kmt),hsj(jmt),hsi(imt)
        do 10 j=1,kmt
          d(j)=1000*j
 10	  continue

	call initiation
	itt=0
	write (6,*) h(21,21)
 900  continue
c marching in time; euler-backward is considered, mxp=0: ordinary l-f;
c mxp=1: only 1st pass of eb is used, mxp=2: complete eb is used

      if(mod(itt,mixeb).ne.0) then
	c2dtts = 2.0d0 * dtts
	c2dtuv = 2.d0 * dtuv
	  do 50 j = 1, jmt
	    do 50 i = 1, imt
	      if(mxp .eq. 0) then
		ub(i,j) = u (i,j)
		vb(i,j) = v (i,j)
		hub(i,j) = hu (i,j)
	      endif
	      u (i,j) = ua(i,j)
	      v (i,j) = va(i,j)
50	      hu (i,j) = hua(i,j)

	    do 60 j = 1, jmt
	      do  60 i = 1, imt
		if(mxp .eq. 0) hb(i,j) = h(i,j)
60		h(i,j) = ha(i,j)
	mxp=0
c        the following is the first pass of euler backward method   
      else
	c2dtts = dtts
	c2dtuv = dtuv
c     if use only first pass of eb, select mxp=1
c     if use both passes of eb, select mxp=2
	mxp = mxp0
	  do 80 j = 1, jmt
	    do 80 i = 1, imt
	      hu (i,j) = hua(i,j)
	      hub(i,j) = hu (i,j)
	      u (i,j) = ua(i,j)
	      ub(i,j) = u (i,j)
	      v (i,j) = va(i,j)
80	      vb(i,j) = v (i,j)

	    do 90 j = 1, jmt
	      do  90 i = 1, imt
		h(i,j) = ha(i,j)
90		hb(i,j) =h(i,j)
      endif

c    after 2nd-pass of euler-backward, come to this pt to continue 
 910  continue
	year=itt*dtts/3600.0d0/24.0d0/365.0d0
	yfact=0.0
 	pi=3.1415926535
	if (year .le. year1) yfact=1.0
        do 915 j=1,jmt
	do 915 i=1,imt
         gp(i,j)=gpri+yfact*gp0(i,j)
         wsx(i,j)=tau(j)+yfact*taud(j)
 915	 continue
      call clinic(itt)
917   call tracer

c        complete the second pass of euler backward method  
      if(mxp.eq.2) then
	mxp = 1
	  do 130 j = 1, jmt
	    do 130 i = 1, imt
	      hu (i,j) = hua(i,j)
	      u (i,j) = ua(i,j)
130	      v (i,j) = va(i,j)
	    do 140 j = 1, jmt
	      do 140 i = 1, imt
140		h (i,j) = ha(i,j)
	go to 910
      endif
	
	if (mod(itt,iwrite) .eq. 0) then
	hmax=0.0d0
	hmin=1.0d+8
	sum=0.0d0
        umax=0
	do 200 j=2,jmm
	do 200 i=2,imm
	sum=sum+h(i,j)-dzi
	if (h(i,j) .gt. hmax) hmax=h(i,j)
	if (h(i,j) .lt. hmin) hmin=h(i,j)
	if (abs(u(i,j)) .gt. umax) umax=abs(u(i,j))
200	continue
c	yfact=1.0
c	if (year .le. p_yr) yfact=year/p_yr
	year=itt*dtts/3600.0d0/24.0d0/365.0d0
	write (6,999)itt*0.001,year,yfact,hmax,hmin,umax,0.01*sum
        write (6,*)  '   '
        hs1=0.0d0
	hs2=0.0d0
	hs3=0.0d0
	do 210 i=2,imm
	do 201 j=2,16
 201	   hs1=hs1+h(i,j)
	do 202 j=2,36
 202	   hs2=hs2+h(i,j)
	do 203 j=2,56
 203	   hs3=hs3+h(i,j)
 210	continue
        tup=13
        tlo=3
        do 220 k=1,kmt
        hct(k)=0.0d0
	do 220 j=2,jmm
	   do 220 i=2,imm
	      hd(k)=h(i,j)/d(k)
	      if (hd(k)>=1) then
                hct(k)=hct(k)+tup*d(k)
              else
                hct(k)=hct(k)+tup*h(i,j)+tlo*(d(k)-h(i,j))
	      end if
 220	      continue
	      do 230 j=2,jmt-1
		 hsj(j)=0
		 do 230 i=2,imt-1
		    hsj(j)=hsj(j)+h(i,j)
 230	continue
 	      do 240 i=2,imt-1
		 hsi(i)=0
		 do 240 j=2,jmt-1
		    hsi(i)=hsi(i)+h(i,j)
 240	continue
       write (201,9999) hs1,hs2,hs3
        write (202,9999) (hct(k),k=1,kmt)
	do 250 j=2,jmm
        write (203,9900) (h(i,j),i=2,imm)
 250	continue
 9900	format (360e20.10)
        write (211,9999) (h(i,31),i=2,imm)
        write (212,9999) (h(i,46),i=2,imm)
        write (213,9999) (h(i,61),i=2,imm)
        write (215,9999) (0.5*(h(i,71)+h(i,72)),i=2,imm)
        write (216,9999) (h(i,82),i=2,imm)
        write (217,9999) (h(i,97),i=2,imm)
        write (218,9999) (h(i,113),i=2,imm)
        write (219,9999) (h(i,128),i=2,imm)
        write (207,9999) (hsj(j),j=2,jmm)
        write (208,9999) (hsi(i),i=2,imm)
        write (200,9999) (yfact)
	end if

      itt   = itt + 1
      if(itt .le. nlast) go to 900
	years=years+nlast*dtts/3600.0d0/24.0d0/365.0d0
      	write (11)  h,u,v
	do 300 j=2,jmm
	do 300 i=2,imm
	write (100,9996) 1.0d-4*h(i,j)
300	continue
	do 310 j=2,jmm1
	yj=float(j-2)/float(jmt-2)
	do 310 i=2,imm1
	xi=float(i-2)/float(imt-2)
	write (101,996) u(i,j)
	write (102,996) v(i,j)
310	continue
	do 320 j=2,jmm
	psij=0.0d0
	do 320 i=imm,2,-1
	psij=psij-h(i,j)*0.25d0*(v(i,j)+v(i-1,j)+
     *                         v(i,j-1)+v(i-1,j-1))*dxt*1.0d-12
	psi(i,j)=psij
320	continue
	do 330 j=2,jmm
	do 330 i=2,imm
330	write (103,996) psi(i,j)

         gravity=980.665
	do 350 j=1,jmt
	do 350 i=1,imt
	eta(i,j)=gpri/gravity*(h(i,j)-dzi)
 350	continue
	do 360 j=2,jmm
	do 360 i=2,imm
360	write (104,996) eta(i,j)
        write (105,997) nlast,iwrite,dzi,dtuv
 997	format (2i10,2f12.2)
998	format (2x,'End of ',f12.2,'years run.')
999	format (10f12.3)
 996	format (10e14.6)
 9996	format (10e20.12)
 9999	format (10e20.10)
      stop 13000
      end

*deck clinic
      subroutine clinic(itt)
	implicit double precision (a-h,o-z)
      parameter(imt=62,jmt=142,imm=imt-1,imm1=imm-1,jmm=jmt-1,
     *     jmm1=jmm-1)
	common /vel/ u(imt,jmt),v(imt,jmt),ua(imt,jmt),va(imt,jmt),
     *  ub(imt,jmt),vb(imt,jmt),h(imt,jmt),ha(imt,jmt),hb(imt,jmt),
     *  hua(imt,jmt),hu(imt,jmt),hub(imt,jmt),f(jmt),gp(imt,jmt),
     *  gp0(imt,jmt),wsx(imt,jmt),tau(jmt),taud(jmt),dh0(imt,jmt)
	common /sca/acor,dxt,dxtr,dtts,dtuv,c2dtts,c2dtuv,drag,gpri,
     *         mxp0,mxp,mixeb,nfirst,nlast,iwrite,dzi,p_yr,year1,year2
	dimension dpdx(imt,jmt),dpdy(imt,jmt),fa(imt,jmt),ga(imt,jmt) 
       do  5 j=2,jmm1
         do  5 i=2,imm1
        dhp=h(i+1,j+1)**2*gp(i,j+1)-h(i,j)**2*gp(i,j)          !g'h^2
        dhq=h(i+1,j)**2*gp(i,j)-h(i,j+1)**2*gp(i,j+1)

        dpdx(i,j)=0.25d0*(dhp+dhq)*dxtr                    !*gpri, 0.5 --> 0.25 
        dpdy(i,j)=0.25d0*(dhp-dhq)*dxtr                    !*gpri
5 	continue
	if (mod(itt,iwrite) .eq. 0) then
        write (6,*) ' '
        end if
c----------compute thicknesses at velocity points
       do 10 j=1,jmm
        do 10 i=1,imm
        hua(i,j)=0.25d0*(h(i,j)+h(i+1,j)+h(i,j+1)+h(i+1,j+1))
10 	continue
c----------add in horizontal viscous effects
      detest=1.d-25
      aauu=0.5d0*am*dxtr*dxtr
        do 20 j=2,jmm1
         do 20 i=2,imm1
       ddip=hu(i+1,j)+hu(i  ,j)
       ddim=hu(i  ,j)+hu(i-1,j)
       ddjp=hu(i,j+1)+hu(i  ,j)
       ddjm=hu(i,j  )+hu(i,j-1)
       fa(i,j)=((u(i+1,j)-u(i  ,j))*ddip
     *           -(u(i  ,j)-u(i-1,j))*ddim
     *           +(u(i,j+1)-u(i  ,j))*ddjp
     *           -(u(i  ,j)-u(i,j-1))*ddjm)*aauu
       ga(i,j)=((v(i+1,j)-v(i  ,j))*ddip
     *           -(v(i  ,j)-v(i-1,j))*ddim
     *           +(v(i,j+1)-v(i  ,j))*ddjp
     *           -(v(i  ,j)-v(i,j-1))*ddjm)*aauu
20 	continue
       do 30 j=2,jmm1
        do 30 i=2,imm1
       fa(i,j)=(fa(i,j)+wsx(i,j)-dpdx(i,j))                 ! delete the hua(i,j) factor
     *    +f(j)*(1.0d0-acor)*hub(i,j)*vb(i,j)
     *    +hub(i,j)*ub(i,j)/c2dtuv
       ga(i,j)=(ga(i,j)	       -dpdy(i,j))
     *    -f(j)*(1.0d0-acor)*hub(i,j)*ub(i,j)
     *    +hub(i,j)*vb(i,j)/c2dtuv
30 	continue
	do 31 j=2,jmm1
	do 31 i=2,imm1
	hk=drag+hub(i,j)/c2dtuv
	hf=hua(i,j)*f(j)*acor
	ddd=1.0d0/(hk**2+hf**2+1.0d-20)
         ua(i,j)=(fa(i,j)*hk+ga(i,j)*hf)*ddd
         va(i,j)=(ga(i,j)*hk-fa(i,j)*hf)*ddd
31 	continue
      	do 40 j=1,jmt
            ua(1,j)=0.0d0
            ua(imm,j)=0.0d0
            ua(imt,j)=0.0d0
            va(1,j)=0.0d0
            va(imm,j)=0.0d0
            va(imt,j)=0.0d0
40 	continue
      	do 50 i=1,imt
            ua(i,1)=0.0d0
            ua(i,jmm)=0.0d0
            ua(i,jmt)=0.0d0
            va(i,1)=0.0d0
            va(i,jmm)=0.0d0
            va(i,jmt)=0.0d0
50 	continue
 1000 return
      end

*deck tracer
      subroutine tracer
	implicit double precision (a-h,o-z)
      parameter(imt=62,jmt=142,imm=imt-1,imm1=imm-1,jmm=jmt-1,
     *     jmm1=jmm-1,ijmt=imt*jmt)
	common /vel/ u(imt,jmt),v(imt,jmt),ua(imt,jmt),va(imt,jmt),
     *  ub(imt,jmt),vb(imt,jmt),h(imt,jmt),ha(imt,jmt),hb(imt,jmt),
     *  hua(imt,jmt),hu(imt,jmt),hub(imt,jmt),f(jmt),gp(imt,jmt),
     *  gp0(imt,jmt),wsx(imt,jmt),tau(jmt),taud(jmt),dh0(imt,jmt)
	common /sca/acor,dxt,dxtr,dtts,dtuv,c2dtts,c2dtuv,drag,gpri,
     *         mxp0,mxp,mixeb,nfirst,nlast,iwrite,dzi,p_yr,year1,year2
      dimension  fux(imt,jmt),fvy(imt,jmt),u0(imt,jmt),v0(imt,jmt)
	dimension had(imt,jmt)
      dimension
     * ppls  (imt,jmt),abf   (imt,jmt),abg  (imt,jmt),
     * af    (imt,jmt),ag    (imt,jmt),
     * dhamn(imt,jmt),dhamx(imt,jmt),
     * pmns  (imt,jmt),qmns  (imt,jmt),rmns (imt,jmt),
     * qpls  (imt,jmt),rpls (imt,jmt),
     * pplsr (imt,jmt),pmnsr (imt,jmt)
	data fux,fvy,u0,v0/ijmt*0.0d0,ijmt*0.0d0,ijmt*0.0d0,ijmt*0.0d0/
	data af,ag/ijmt*0.0d0,ijmt*0.0d0/
       do 10 j = 2, jmm
	do 10 i = 2,imm1
	   u0(i,j)=0.5d0*(ua(i,j) + ua(i,j-1))
	fux(i,j)=0.5d0*((h(i,j)+h(i+1,j))*u0(i,j)+
     *               (h(i,j)-h(i+1,j))*dabs(u0(i,j)))
10	continue

      do 20 i = 2, imm
	do 20 j = 2,jmm1
	  v0(i,j)= 0.5d0*(va(i,j)+va(i-1,j))
	fvy(i,j)=0.5d0*((h(i,j)+h(i,j+1))*v0(i,j)+
     *		     (h(i,j)-h(i,j+1))*dabs(v0(i,j)))
20	continue

      fact=dxtr*c2dtts
      do 30 j=2,jmm
         do 30 i=2,imm
            had(i,j)=hb(i,j)-
     *      (fux(i,j)-fux(i-1,j)+fvy(i,j)-fvy(i,j-1))*fact
30 	continue
      do 40 j=1,jmt
            had(1,j)=had(2,j)
40            had(imt,j)=had(imm,j)
      do 50 i=1,imt
            had(i,1)=had(i,2)
50            had(i,jmt)=had(i,jmm)
c
c----------compute 'anti-diffusive' fluxes (high-order minus low-order)
c
       do 130 j=2,jmm
         do 130 i=2,imm1
          af(i,j)=
     *  (0.5d0*(u(i,j-1)+u(i,j))*(.5625d0*(h(i+1,j)+h(i,j))
     *       -.0625d0*(h(i+2,j)+h(i-1,j)))-fux(i,j))
130	continue
	do 135 i=2,imm
	do 135 j=2,jmm1
          ag(i,j)=
     *  (0.5d0*(v(i,j)+v(i-1,j))*(.5625d0*(h(i,j+1)+h(i,j))
     *       -.0625d0*(h(i,j+2)+h(i,j-1)))-fvy(i,j))
135 	continue
cc
       do 140 j=1,jmt
        do 140 i=1,imt
         abf(i,j)=dabs(af(i,j))
         abg(i,j)=dabs(ag(i,j))
140 	continue
c
c----------compute maximum positive (negative) change in dzt
c
       do 170 j=2,jmm
        do 170 i=2,imm
        dp=had(i-1,j)+had(i,j+1)
        dm=dabs(had(i-1,j)-had(i,j+1))
        dq=had(i,j-1)+had(i+1,j)
        dn=dabs(had(i,j-1)-had(i+1,j))
         temp1=dp+dm
         temp2=dq+dn
         temp3=0.25d0*(temp1+temp2+dabs(temp1-temp2))
         dhamx(i,j)=0.5d0*(had(i,j)+temp3+dabs(had(i,j)-temp3))
         temp4=dp-dm
         temp5=dq-dn
         temp6=0.25d0*(temp4+temp5-dabs(temp5-temp4))
         dhamn(i,j)=0.5d0*(had(i,j)+temp6-dabs(had(i,j)-temp6))
170 	continue
c
c----------compute ratio of maximum positive (negative) change in dzt
c           to sum of fluxes into (out of) grid box
c
       do 80 j=2,jmm
         do 80 i=2,imm
            temp1=af(i-1,j)+abf(i-1,j)
            temp2=af(i,j)-abf(i,j)
            temp3=ag(i,j-1)+abg(i,j-1)
            temp4=ag(i,j)-abg(i,j)
            ppls(i,j)=temp1-temp2+temp3-temp4
            temp5=af(i,j)+abf(i,j)
            temp6=af(i-1,j)-abf(i-1,j)
            temp7=ag(i,j)+abg(i,j)
            temp8=ag(i,j-1)-abg(i,j-1)
            pmns(i,j)=temp5-temp6+temp7-temp8
            qpls(i,j)=dhamx(i,j)-had (i,j)
            qmns(i,j)=had (i,j)-dhamn(i,j)
   80 continue
      do 85 j=2,jmm
      do 85 i=2,imm
       pplsr(i,j)=0.0d0
       if (ppls(i,j) .ne. 0.0) pplsr(i,j)=1.0d0/ppls(i,j)
   85 continue
      do 90 j=2,jmm 

      do 90 i=2,imm
       pmnsr(i,j)=0.0d0
       if (pmns(i,j) .ne. 0.0) pmnsr(i,j)=1.0d0/pmns(i,j)
   90 continue
         fact=dxt/c2dtts
       do 100 j=2,jmm
         do 100 i=2,imm
            temp1=qpls(i,j)*pplsr(i,j)*fact
            temp2=qmns(i,j)*pmnsr(i,j)*fact
            rpls(i,j)=0.5d0+temp1-dabs(0.5d0-temp1)
            rmns(i,j)=0.5d0+temp2-dabs(0.5d0-temp2)
  100 continue

       do 110 j=1,jmt
            rpls(1,j)=rpls(2,j)
            rmns(1,j)=rmns(2,j)
            rpls(imt,j)=rpls(imm,j)
            rmns(imt,j)=rmns(imm,j)
  110 continue
       do 120 i=1,imt
            rpls(i,1)=rpls(i,2)
            rmns(i,1)=rmns(i,2)
            rpls(i,jmt)=rpls(i,jmm)
            rmns(i,jmt)=rmns(i,jmm)
  120 continue
c
c----------compute coefficients that determine contributions from
c           low-order and high-order fluxes
c
       do 220 j=2,jmm
        do 220 i=2,imm
         temp1=rpls(i+1,j)+rmns(i,j)-dabs(rmns(i,j)-rpls(i+1,j))
         temp2=rpls(i,j)+rmns(i+1,j)-dabs(rmns(i+1,j)-rpls(i,j))
       af(i,j)=.25d0*((temp1+temp2)*af(i,j)+(temp1-temp2)*abf(i,j))
         temp3=rpls(i,j+1)+rmns(i,j)-dabs(rmns(i,j)-rpls(i,j+1))
         temp4=rpls(i,j)+rmns(i,j+1)-dabs(rmns(i,j+1)-rpls(i,j))
       ag(i,j)=.25d0*((temp3+temp4)*ag(i,j)+(temp3-temp4)*abg(i,j))
  220 continue
c
c----------apply correction to dzta
c
      fact=dxtr*c2dtts
       do 260 j=2,jmm
         do 260 i=2,imm
        ha(i,j)=(had(i,j)-(af(i,j)-af(i-1,j)
     *       +ag(i,j)-ag(i,j-1))*fact)
        ha(i,j)=max(ha(i,j),0.0d0)
  260 continue

	sum1=0.0d0
	do 2000 j=2,jmm
	do 2000 i=2,imm
	sum1=sum1+ha(i,j)
2000	continue
cc	write (6,998) year,sum*1.0d-8,sum1*1.0d-8
998 	format (3f12.6)

      return
      end

*deck initiation
	subroutine initiation
	implicit double precision (a-h,o-z)
      parameter(imt=62,jmt=142,imm=imt-1,imm1=imm-1,jmm=jmt-1,
     *     jmm1=jmm-1)
	common /vel/ u(imt,jmt),v(imt,jmt),ua(imt,jmt),va(imt,jmt),
     *  ub(imt,jmt),vb(imt,jmt),h(imt,jmt),ha(imt,jmt),hb(imt,jmt),
     *  hua(imt,jmt),hu(imt,jmt),hub(imt,jmt),f(jmt),gp(imt,jmt),
     *  gp0(imt,jmt),wsx(imt,jmt),tau(jmt),taud(jmt),dh0(imt,jmt)
	common /sca/acor,dxt,dxtr,dtts,dtuv,c2dtts,c2dtuv,drag,gpri,
     *         mxp0,mxp,mixeb,nfirst,nlast,iwrite,dzi,p_yr,year1,year2
        dimension ylat(jmt),yarc(jmt),fh(jmt)
	open (15,file='input_dv2')
      read  (15,*) nfirst,nlast,iwrite,am,dzi,dtts,dtuv,dxt
      write (6,61) nfirst,nlast,iwrite,am,0.01*dzi,dtts,dtuv,dxt
      read  (15,*) drag,gpri,tau_fact,h_scale
      write (6,62) drag,gpri,tau_fact,h_scale
      read  (15,*)  phi0,dphi,d_tau,dgpri,year1,year2,ddh,iwe,iea
      write (6,62)  phi0,dphi,d_tau,dgpri,year1,year2,ddh,iwe,iea
 61	format (3i8,e12.2,3f12.0,e12.2 )
 62	format (7f10.2,3i8)
	open (112,file='profile.data')
	open (113,file='dgp.data')
c----------define physical constants
      pi=3.1415926
      pi8=pi/180
      acor=0.6666667d0
      grav=980.6d0
      omega=3.14159265d0/43082.0d0
      fzero=8.3652d-5       !central latitude at 35N
      beta=2.2367d-13
      mxp0 = 2
      mxp  = 0
      mixeb=23
c----------compute coriolis parameter (central latitude = 45 degrees)
      jfact=(jmt-2)/2+35
       do 5 j=1,jmt
            fh(j)=fzero+beta*(float(j-jfact))*dxt
            f(j) =fzero+beta*(float(j-jfact)-0.5)*dxt
	write (111,99) f(j)*1.0e5,fh(j)*1.0e5,float(j)
99	format (5f12.6)
5 	continue
      dzir=1.0d0/dzi
      dxtr=1.0d0/dxt
      dxt2r=.5d0*dxtr
      gpr1=gpri*dxt2r
c----------compute wind stress factor
	do 6 j=1,jmt
        do 6 i=1,imt
           dh0(i,j)=0
 6	   gp0(i,j)=0
c             js=81
c	     jn=110
              js=111
              jn=141
c              is=23
c              in=65
	do 7 j=js,jn
        do 7 i=iwe,iea
              ddphj=pi*float(j-js)/float(jn-js)
              ddphi=pi*float(i-iwe)/float(iea-iwe)
 7	   dh0(i,j)=ddh*sin(ddphi)*sin(ddphj)
	do 8 j=1,jmt
        do 8 i=1,imt
 8	  write (113,993) dh0(i,j)
 993	  format (f12.3)

      do 10 j=1,jmt
       ylat(j)=float(j)-71
       yarcj=(abs(float(j)-71))*pi8
       tau(j)=0.2-0.8*sin(6.0d0*yarcj)
     *        -0.5d0*(1.0d0-tanh(10.0d0*yarcj))
     *        -0.5d0*(1.0d0-tanh(10.0d0*(0.5*pi-yarcj)))
	yyy=(ylat(j)-phi0)/dphi
	taud(j)=d_tau*exp(-yyy**2)
 10	continue

         do 11 j=1,jmt
          do 11 i=1,imt
            wsx(i,j)=tau(j)+taud(j)
 11	 continue
         do 12 j=1,jmt
            write (112,99) float(j),f(j),wsx(1,j),tau(j),gpri+gp0(20,j)
12 	continue
c     end of the note
	if (nfirst .eq. 0) then
	do 20 j=1,jmt
	do 20 i=1,imt
	ha(i,j)=dzi
	hua(i,j)=dzi
	h(i,j)=dzi
	hu(i,j)=dzi
	hb(i,j)=dzi
	hub(i,j)=dzi
	u(i,j)=0.0d0
	v(i,j)=0.0d0
	ub(i,j)=0.0d0
	vb(i,j)=0.0d0
	ua(i,j)=0.0d0
	va(i,j)=0.0d0
20	continue
	else
	read (12) h,u,v
	do 120 j=1,jmt
	do 120 i=1,imt
	ha (i,j)=h(i,j)+dh0(i,j)
	hua(i,j)=h(i,j)+dh0(i,j)
	hu (i,j)=h(i,j)+dh0(i,j)
	hb (i,j)=h(i,j)+dh0(i,j)
	hub(i,j)=h(i,j)+dh0(i,j)
	ub(i,j)=u(i,j)
	vb(i,j)=v(i,j)
	ua(i,j)=u(i,j)
	va(i,j)=v(i,j)
120	continue
	write (6,*) h(15,125),dh0(15,125),ha(15,125)
	end if

	return
	end

