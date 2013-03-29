c
c*********************************************************************
c*********************************************************************
c***                                                               ***
c***    " NIST Extension Temperature and Density Multiphase "      ***
c***                                                               ***
c***                   VERSION  # 1.0.0 02/05/13                   ***
c***                                                               ***
c***        Copyright 2013 CSEWorks INC All Rights Reserved        ***
c***                   PROPRIETARY INFORMATION                     ***
c***                For CSEWorks INC Use Only.                     ***
c***                                                               ***
c*********************************************************************
c*********************************************************************
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDETDFLSH                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDETDFLSH: Extended TDFLSH
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDETDFLSH (t,D,z,p,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDETDFLSH
c     dll_export TDETDFLSH
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        p=D*R*t
        Dl=D
        Dv=D
        call THERM (t,D,z,ptherm,e,h,s,cv,cp,w,hjt)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDETDFLSH warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDETDFLSH error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        q=999.d0     !quality undefined
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDETDFLSH error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call THERM (t,D,z,p,e,h,s,cv,cp,w,hjt)
        Dl=D
        Dv=D
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call THERM (t,D,z,p,e,h,s,cv,cp,w,hjt)
        Dl=D
        Dv=D
        if (p.le.pc) then
          q=998.d0     !superheated, but cannot define quality
        else
          q=999.d0     !quality undefined
        end if
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDETDFLSH error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          D=0.0d0
          Dl=0.0d0
          Dv=0.0d0
          q=999.d0     !quality undefined
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[ETPFLSH error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            D=0.0d0
            Dl=0.0d0
            Dv=0.0d0
            q=999.d0     !quality undefined
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call THERM (t,D,z,p,e,h,s,cv,cp,w,hjt)
          Dl=D
          Dv=D
          call ENTHAL (t,Dvdew,z,hvdew)
c         call ENTHAL (t,Dlbub,z,hlbub)
c         q=(h-hlbub)/(hvdew-hlbub)
c  compute quality based on volumes (possible for T,h to be
c  double-valued in compressed liquid)
          if (abs(dvdew-dlbub).gt.1.d-12)
     &        q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDETDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call THERM (t,Dl,x,ptherm,el,hl,sl,cvl,cpl,wl,hjt)
          call THERM (t,Dv,y,ptherm,ev,hv,sv,cvv,cpv,wv,hjt)
          if (nc.ge.2) p=ptherm
          e=q*ev+(1.0d0-q)*el
          h=q*hv+(1.0d0-q)*hl
          s=q*sv+(1.0d0-q)*sl
c  compute the harmonical average for 2-phase states
          tmp=(q*Dv+(1.0d0-q)*Dl)
          w=tmp*(q/(Dv*wv*wv) + (1.0d0-q)/(Dl*wl*wl))
          w=SQRT(1.0d0/w);
          cp=tmp*(q/(Dv*cpv*cpv) + (1.0d0-q)/(Dl*cpl*cpl))
          cp=SQRT(1.0d0/cp);
          cv=tmp*(q/(Dv*cvv*cvv) + (1.0d0-q)/(Dl*cvl*cvl))
          cv=SQRT(1.0d0/cv);
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDETDFLSH
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDETDFLSH2                        ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDETDFLSH2: Extended TDFLSH
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDETDFLSH2 (t,D,z,x,y,thermv,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDETDFLSH2
c     dll_export TDETDFLSH2
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension thermv(29)
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      do i=1,29
        thermv(i) = 0.d0
      enddo
      ierr = 0
      herr = ' '
      p    = 0.d0
      Dl   = 0.d0
      Dv   = 0.d0
      q    = 0.d0
      e    = 0.d0
      h    = 0.d0
      s    = 0.d0
      cv   = 0.d0
      cp   = 0.d0
      w    = 0.d0
      zc   = 0.d0
      hjt  = 0.d0
      a    = 0.d0
      g    = 0.d0
      xkappa   = 0.d0
      beta     = 0.d0
      dpdrho   = 0.d0
      dpt      = 0.d0
      drhodt   = 0.d0
      drhodp   = 0.d0
      d2pdd2   = 0.d0
      d2pdt2   = 0.d0
      d2pdtrho = 0.d0
      dhdt_d   = 0.d0
      dhdt_p   = 0.d0
      dhdd_t   = 0.d0
      dhdd_p   = 0.d0
      dhdp_t   = 0.d0
      dhdp_d   = 0.d0
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        p=D*R*t
        Dl=D
        Dv=D
        call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
c  Finally return the thermodynamic properties
        thermv(1)  = p
        thermv(2)  = Dl
        thermv(3)  = Dv
        thermv(4)  = q
        thermv(5)  = e
        thermv(6)  = h
        thermv(7)  = s
        thermv(8)  = cv
        thermv(9)  = cp
        thermv(10) = w
        thermv(11) = zc
        thermv(12) = hjt
        thermv(13) = a
        thermv(14) = g
        thermv(15) = xkappa
        thermv(16) = beta
        thermv(17) = dpdrho
        thermv(18) = dpt
        thermv(19) = drhodt
        thermv(20) = drhodp
        thermv(21) = d2pdd2
        thermv(22) = d2pdt2
        thermv(23) = d2pdtrho
        thermv(24) = dhdt_d
        thermv(25) = dhdt_p
        thermv(26) = dhdd_t
        thermv(27) = dhdd_p
        thermv(28) = dhdp_t
        thermv(29) = dhdp_d
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDETDFLSH2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDETDFLSH2 error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        q=999.d0     !quality undefined
        ierr=ierrl
        thermv(4) = q
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDETDFLSH2 error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
        thermv(4) = q
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then                 !if 0
c  temperature less than triple point temperature, but in the gas phase.
        call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        Dl=D
        Dv=D
      elseif (t.ge.tc .or. D.lt.1.0d-10) then                !elseif 0
c  super-critical state or rho = 0 (x = y = z as set above)
        call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        Dl=D
        Dv=D
        if (p.le.pc) then
          q=998.d0     !superheated, but cannot define quality
        else
          q=999.d0     !quality undefined
        end if
      else                                                     !else 0
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDETDFLSH2 error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          D=0.0d0
          Dl=0.0d0
          Dv=0.0d0
          q=999.d0     !quality undefined
          thermv(2) = Dl
          thermv(3) = Dv
          thermv(4) = q
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[ETPFLSH error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            D=0.0d0
            Dl=0.0d0
            Dv=0.0d0
            q=999.d0     !quality undefined
            thermv(2) = Dl
            thermv(3) = Dv
            thermv(4) = q
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then                     !if 1
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          Dl=D
          Dv=D
          call ENTHAL (t,Dvdew,z,hvdew)
c         call ENTHAL (t,Dlbub,z,hlbub)
c         q=(h-hlbub)/(hvdew-hlbub)
c  compute quality based on volumes (possible for T,h to be
c  double-valued in compressed liquid)
          if (abs(dvdew-dlbub).gt.1.d-12)
     &        q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
        else                                                   !else 1
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              thermv(4) = p
              thermv(2) = Dl
              thermv(3) = Dv
              thermv(4) = q
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDETDFLSH2--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call THERM2 (t,Dl,x,ptherm,el,hl,sl,cvl,cpl,wl,zcl,hjtl,al,gl,
     &      xkappal,betal,dpdrhol,d2pdd2l,dptl,drhodtl,drhodpl,d2pdt2l,
     &      d2pdtrhol,spare3,spare4)
c
          call THERM2 (t,Dv,y,ptherm,ev,hv,sv,cvv,cpv,wv,zcv,hjtv,av,gv,
     &      xkappav,betav,dpdrhov,d2pdd2v,dptv,drhodtv,drhodpv,d2pdt2v,
     &      d2pdtrhov,spare3,spare4)
c
          call DHD1 (t,Dl,x,dhdt_dl,dhdt_pl,dhdd_tl,dhdd_pl,
     &      dhdp_tl,dhdp_dl)
c
          call DHD1 (t,Dv,y,dhdt_dv,dhdt_pv,dhdd_tv,dhdd_pv,
     &      dhdp_tv,dhdp_dv)
c          
          if (nc.ge.2) p=ptherm
c  bulk properties are weighted average of liquid and vapor phases
          alpha=1.0d0-q               !alpha is liq fraction,
          e=alpha*el+q*ev
          h=alpha*hl+q*hv
          s=alpha*sl+q*sv
          zc=alpha*zcl+q*zcv
          hjt=alpha*hjtl+q*hjtv
          a=alpha*al+q*av
          g=alpha*gl+q*gv
          xkappa=alpha*xkappal+q*xkappav
          beta=alpha*betal+q*betav
          dpdrho=alpha*dpdrhol+q*dpdrhov
          d2pdd2=alpha*d2pdd2l+q*d2pdd2v
          drhodt=alpha*drhodtl+q*drhodtv
          dpt=-dpdrho*drhodt
          drhodp=1.0d0/dpdrho
          d2pdt2=alpha*d2pdt2l+q*d2pdt2v
          d2pdtrho=alpha*d2pdtrhol+q*d2pdtrhov
c
          dhdt_d=alpha*dhdt_dl+q*dhdt_dv
          dhdt_p=alpha*dhdt_pl+q*dhdt_pv
          dhdd_t=alpha*dhdd_tl+q*dhdd_tv
          dhdd_p=alpha*dhdd_pl+q*dhdd_pv
          dhdp_t=alpha*dhdp_tl+q*dhdp_tv
          dhdp_d=alpha*dhdp_dl+q*dhdp_dv
c  compute the harmonical average for 2-phase states
          tmp=(q*Dv+alpha*Dl)
          w=tmp*(q/(Dv*wv*wv) + alpha/(Dl*wl*wl))
          w=SQRT(1.0d0/w);
          cp=tmp*(q/(Dv*cpv*cpv) + alpha/(Dl*cpl*cpl))
          cp=SQRT(1.0d0/cp);
          cv=tmp*(q/(Dv*cvv*cvv) + alpha/(Dl*cvl*cvl))
          cv=SQRT(1.0d0/cv);          
        end if                                               !end if 1
      end if                                                 !end if 0
c
c  Finally return the thermodynamic properties
      thermv(1)  = p
      thermv(2)  = Dl
      thermv(3)  = Dv
      thermv(4)  = q
      thermv(5)  = e
      thermv(6)  = h
      thermv(7)  = s
      thermv(8)  = cv
      thermv(9)  = cp
      thermv(10) = w
      thermv(11) = zc
      thermv(12) = hjt
      thermv(13) = a
      thermv(14) = g
      thermv(15) = xkappa
      thermv(16) = beta
      thermv(17) = dpdrho
      thermv(18) = dpt
      thermv(19) = drhodt
      thermv(20) = drhodp
      thermv(21) = d2pdd2
      thermv(22) = d2pdt2
      thermv(23) = d2pdtrho
      thermv(24) = dhdt_d
      thermv(25) = dhdt_p
      thermv(26) = dhdd_t
      thermv(27) = dhdd_p
      thermv(28) = dhdp_t
      thermv(29) = dhdp_d
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                       !subroutine TDETDFLSH2
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDESSOUND                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDESSOUND: Extended SSOUND
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDESSOUND (t,D,z,w,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDESSOUND
c     dll_export TDESSOUND
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call THERM (t,D,z,ptherm,e,h,s,cv,cp,w,hjt)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDESSOUND warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDESSOUND error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDESSOUND error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call THERM (t,D,z,p,e,h,s,cv,cp,w,hjt)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call THERM (t,D,z,p,e,h,s,cv,cp,w,hjt)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDESSOUND error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDESSOUND error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call THERM (t,D,z,p,e,h,s,cv,cp,w,hjt)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDESSOUND--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call THERM (t,Dl,x,pl,el,hl,sl,cvl,cpl,wl,hjt)
          call THERM (t,Dv,y,pv,ev,hv,sv,cvv,cpv,wv,hjt)
          if (nc.ge.2) p=0.5d0*(pl+pv)
c  compute the harmonical average for 2-phase states
          tmp=(q*Dv+(1.0d0-q)*Dl)
          w=tmp*(q/(Dv*wv*wv) + (1.0d0-q)/(Dl*wl*wl))
          w=SQRT(1.0d0/w);
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDESSOUND      
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEPRESS                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEPRESS: Extended PRESS
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEPRESS (t,D,z,p,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEPRESS
c     dll_export TDEPRESS
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        p=D*R*t
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEPRESS warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEPRESS error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEPRESS error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEPRESS error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEPRESS error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEPRESS--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                         !subroutine TDEPRESS
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEENTRO                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEENTRO: Extended ENTRO
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEENTRO (t,D,z,s,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEENTRO
c     dll_export TDEENTRO
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call ENTRO (t,D,z,s)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEENTRO warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEENTRO error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEENTRO error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call ENTRO (t,D,z,s)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call ENTRO (t,D,z,s)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEENTRO error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEENTRO error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call ENTRO (t,D,z,s)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEENTRO--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call ENTRO (t,Dl,x,sl)
          call ENTRO (t,Dv,y,sv)
          s=q*sv+(1.0d0-q)*sl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                         !subroutine TDEENTRO
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEENTHAL                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEENTHAL: Extended ENTHAL
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEENTHAL (t,D,z,h,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEENTHAL
c     dll_export TDEENTHAL
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call ENTHAL (t,D,z,h)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEENTHAL warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEENTHAL error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEENTHAL error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call ENTHAL (t,D,z,h)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call ENTHAL (t,D,z,h)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEENTHAL error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEENTHAL error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call ENTHAL (t,D,z,h)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEENTHAL--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call ENTHAL (t,Dl,x,hl)
          call ENTHAL (t,Dv,y,hv)
          h=q*hv+(1.0d0-q)*hl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDEENTHAL
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEENERGY                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEENERGY: Extended ENERGY
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEENERGY (t,D,z,e,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEENERGY
c     dll_export TDEENERGY
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call ENERGY (t,D,z,e)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEENERGY warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEENERGY error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEENERGY error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call ENERGY (t,D,z,e)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call ENERGY (t,D,z,e)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEENERGY error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEENERGY error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call ENERGY (t,D,z,e)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEENERGY--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call ENERGY (t,Dl,x,el)
          call ENERGY (t,Dv,y,ev)
          e=q*ev+(1.0d0-q)*el
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDEENERGY
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDECVCP                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDECVCP: Extended CVCP
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDECVCP (t,D,z,cv,cp,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDECVCP
c     dll_export TDECVCP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call CVCP (t,D,z,cv,cp)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDECVCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDECVCP error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDECVCP error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call CVCP (t,D,z,cv,cp)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call CVCP (t,D,z,cv,cp)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDECVCP error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDECVCP error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call CVCP (t,D,z,cv,cp)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDECVCP--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call CVCP (t,Dl,x,cvl, cpl)
          call CVCP (t,Dv,y,cvv, cpv)
c  compute the harmonical average for 2-phase states
          tmp=(q*Dv+(1.0d0-q)*Dl)
          cp=tmp*(q/(Dv*cpv*cpv) + (1.0d0-q)/(Dl*cpl*cpl))
          cp=SQRT(1.0d0/cp);
          cv=tmp*(q/(Dv*cvv*cvv) + (1.0d0-q)/(Dl*cvl*cvl))
          cv=SQRT(1.0d0/cv);
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                          !subroutine TDECVCP
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDPDD                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDPDD: Extended DPDD
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDPDD (t,D,z,dpdrho,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDPDD
c     dll_export TDEDPDD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DPDD (t,D,z,dpdrho)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDPDD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDPDD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDPDD error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DPDD (t,D,z,dpdrho)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DPDD (t,D,z,dpdrho)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDPDD error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDPDD error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DPDD (t,D,z,dpdrho)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDPDD--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DPDD (t,Dl,x,dpdrhol)
          call DPDD (t,Dv,y,dpdrhov)
          dpdrho=q*dpdrhov+(1.0d0-q)*dpdrhol
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                          !subroutine TDEDPDD
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDPDT                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDPDT: Extended DPDT
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDPDT (t,D,z,dpt,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDPDT
c     dll_export TDEDPDT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DPDT (t,D,z,dpt)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDPDT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDPDT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDPDT error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DPDT (t,D,z,dpt)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DPDT (t,D,z,dpt)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDPDT error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDPDT error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DPDT (t,D,z,dpt)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDPDT--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DPDT (t,Dl,x,dptl)
          call DPDT (t,Dv,y,dptv)
          dpt=q*dptv+(1.0d0-q)*dptl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                          !subroutine TDEDPDT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDDDT                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDDDT: Extended DDDT
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDDDT (t,D,z,drhodt,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDDDT
c     dll_export TDEDDDT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DDDT (t,D,z,drhodt)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDDDT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDDDT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDDDT error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DDDT (t,D,z,drhodt)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DDDT (t,D,z,drhodt)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDDDT error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDDDT error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DDDT (t,D,z,drhodt)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDDDT--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DDDT (t,Dl,x,drhodtl)
          call DDDT (t,Dv,y,drhodtv)
          drhodt=q*drhodtv+(1.0d0-q)*drhodtl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                          !subroutine TDEDDDT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDPDD2                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDPDD2: Extended DPDD2
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDPDD2 (t,D,z,d2pdrho2,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDPDD2
c     dll_export TDEDPDD2
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DPDD2 (t,D,z,d2pdrho2)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDPDD2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDPDD2 error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDPDD2 error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DPDD2 (t,D,z,d2pdrho2)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DPDD2 (t,D,z,d2pdrho2)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDPDD2 error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDPDD2 error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DPDD2 (t,D,z,d2pdrho2)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDPDD2--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DPDD2 (t,Dl,x,d2pdrho2l)
          call DPDD2 (t,Dv,y,d2pdrho2v)
          d2pdrho2=q*d2pdrho2v+(1.0d0-q)*d2pdrho2l
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                         !subroutine TDEDPDD2
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDPDT2                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDPDT2: Extended DPDT2
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDPDT2 (t,D,z,d2pdt2,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDPDT2
c     dll_export TDEDPDT2
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDPDT2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDPDT2 error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDPDT2 error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDPDT2 error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDPDT2 error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDPDT2--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call THERM2 (t,Dl,z,pl,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2l,d2pdtrho,
     &      spare3,spare4)
          call THERM2 (t,Dv,z,pv,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2v,d2pdtrho,
     &      spare3,spare4)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          d2pdt2=q*d2pdt2v+(1.0d0-q)*d2pdt2l
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                         !subroutine TDEDPDT2
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDPDTD                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDPDTD: Extended DPDTD
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDPDTD (t,D,z,d2pdtrho,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDPDTD
c     dll_export TDEDPDTD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDPDTD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDPDTD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDPDTD error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDPDTD error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDPDTD error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call THERM2 (t,D,z,p,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDPDTD--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call THERM2 (t,Dl,z,pl,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrhol,
     &      spare3,spare4)
          call THERM2 (t,Dv,z,pv,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrhov,
     &      spare3,spare4)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          d2pdtrho=q*d2pdtrhov+(1.0d0-q)*d2pdtrhol
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                         !subroutine TDEDPDTD
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDHDTCD                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDHDTCD: Extended DHDTCD
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDHDTCD (t,D,z,dhdt_d,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDHDTCD
c     dll_export TDEDHDTCD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDHDTCD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDHDTCD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDHDTCD error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDHDTCD error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDHDTCD error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDHDTCD--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DHD1 (t,Dl,x,dhdt_dl,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          call DHD1 (t,Dv,y,dhdt_dv,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          dhdt_d=q*dhdt_dv+(1.0d0-q)*dhdt_dl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDEDHDTCD
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDHDTCP                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDHDTCP: Extended DHDTCP
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDHDTCP (t,D,z,dhdt_p,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDHDTCP
c     dll_export TDEDHDTCP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDHDTCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDHDTCP error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDHDTCP error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDHDTCP error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDHDTCP error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDHDTCP--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DHD1 (t,Dl,x,dhdt_d,dhdt_pl,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          call DHD1 (t,Dv,y,dhdt_d,dhdt_pv,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          dhdt_p=q*dhdt_pv+(1.0d0-q)*dhdt_pl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDEDHDTCP
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDHDDCT                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDHDDCT: Extended DHDDCT
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDHDDCT (t,D,z,dhdd_t,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDHDDCT
c     dll_export TDEDHDDCT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDHDDCT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDHDDCT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDHDDCT error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDHDDCT error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDHDDCT error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDHDDCT--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DHD1 (t,Dl,x,dhdt_d,dhdt_p,dhdd_tl,dhdd_p,dhdp_t,dhdp_d)
          call DHD1 (t,Dv,y,dhdt_d,dhdt_p,dhdd_tv,dhdd_p,dhdp_t,dhdp_d)
          dhdd_t=q*dhdd_tv+(1.0d0-q)*dhdd_tl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDEDHDDCT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDHDDCP                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDHDDCP: Extended DHDDCP
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDHDDCP (t,D,z,dhdd_p,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDHDDCP
c     dll_export TDEDHDDCP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDHDDCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDHDDCP error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDHDDCP error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDHDDCP error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDHDDCP error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDHDDCP--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DHD1 (t,Dl,x,dhdt_d,dhdt_p,dhdd_t,dhdd_pl,dhdp_t,dhdp_d)
          call DHD1 (t,Dv,y,dhdt_d,dhdt_p,dhdd_t,dhdd_pv,dhdp_t,dhdp_d)
          dhdd_p=q*dhdd_pv+(1.0d0-q)*dhdd_pl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDEDHDDCP
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDHDPCT                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDHDPCT: Extended DHDPCT
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDHDPCT (t,D,z,dhdp_t,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDHDPCT
c     dll_export TDEDHDPCT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDHDPCT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDHDPCT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDHDPCT error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDHDPCT error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDHDPCT error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDHDPCT--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DHD1 (t,Dl,x,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_tl,dhdp_d)
          call DHD1 (t,Dv,y,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_tv,dhdp_d)
          dhdp_t=q*dhdp_tv+(1.0d0-q)*dhdp_tl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDEDHDPCT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDEDHDPCD                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDEDHDPCD: Extended DHDPCD
c  
c  Extended flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDEDHDPCD (t,D,z,dhdp_d,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDEDHDPCD
c     dll_export TDEDHDPCD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr1,herr2,herrl
      dimension x(ncmax),y(ncmax),z(ncmax)
      dimension xbub(ncmax),xdew(ncmax),ybub(ncmax),ydew(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
c
c  set output liquid and vapor compositions to input values
c  zero output values for undefined components
      if (icomp.eq.0) then
        do i=1,nc
          x(i)=z(i)
          y(i)=z(i)
        enddo
      else
        x(icomp)=1.d0
        y(icomp)=1.d0
      endif
c     if (nc.lt.ncmax) then
c       do i=nc+1,ncmax
c         x(i)=0.0d0
c         y(i)=0.0d0
c       enddo
c     end if
      if (d.lt.1.0d-14) then
c  ideal-gas
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        RETURN
      endif
c
c  check that input conditions are within limits
c
      pdum=0.0d0
      call LIMITX ('EOS',t,D,pdum,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[TDEDHDPCD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDEDHDPCD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDEDHDPCD error 220] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      if (icomp.ne.0 .and. t.lt.ttp(icomp)) then
c  temperature less than triple point temperature, but in the gas phase.
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      elseif (t.ge.tc .or. D.lt.1.0d-10) then
c  super-critical state or rho = 0 (x = y = z as set above)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        call PRESS (t,D,z,p)
c  Dl and Dv not needed
      else
c  sub-critical state--call saturation routine to determine phase
        call SATT (t,z,2,pdew,Dldew,Dvdew,xdew,ydew,ierr,herr1)
        if (ierr.ne.0) then
          ierr=222
          write (herr,1222) herr1(1:193),hnull
 1222     format ('[TDEDHDPCD error 222] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        if (icomp.ne.0 .and. ianc(icomp).eq.0) then
c  if pure, bubble point density provided by call to SATT at dew point
          Dlbub=Dldew
          Dvbub=Dvdew
          pbub=pdew
        else
c  if mixture, call saturation routine again at bubble point
          call SATT (t,z,1,pbub,Dlbub,Dvbub,xbub,ybub,ierr,herr1)
          if (ierr.ne.0) then
            ierr=221
            write (herr,1221) herr1(1:190),hnull
 1221       format ('[TDEDHDPCD error 221] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
        end if
        if (D.le.Dvdew .or. D.ge.Dlbub) then
c  single-phase (liq or vapor) (pure or mixture) (x = y = z as set above)
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          call PRESS (t,D,z,p)
c  Dl and Dv not needed
        else
          if (icomp.ne.0) then
c  special case:  pure-fluid two-phase (x = y = z as set above)
            if (abs(dvdew-dlbub).gt.1.d-12)
     &          q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
            p=q*pdew+(1.0d0-q)*pbub
            Dl=Dlbub
            Dv=Dvdew
          else
c  general case:  mixture 2-phase
            ksat=1   !bubble and dew point data provided to TDFL2
            call TDFL2 (t,D,z,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &                  p,Dl,Dv,x,y,q,ierr,herr)
            if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by TDFL2
              q=999.d0     !quality undefined
              RETURN
            end if
c           write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225       format (1x,' TDEDHDPCD--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          call DHD1 (t,Dl,x,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_dl)
          call DHD1 (t,Dv,y,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_dv)
          dhdp_d=q*dhdp_dv+(1.0d0-q)*dhdp_dl
        end if
      end if
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierr,herrl)
      if (ierr.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierr,herrl(1:234),hnull
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierr.gt.0) then
        write (herr,1007) ierr,herrl(1:236),hnull
        call ERRMSG (ierr,herr)
        q=999.d0     !quality undefined
      end if
c
      RETURN
      end                                        !subroutine TDEDHDPCD

