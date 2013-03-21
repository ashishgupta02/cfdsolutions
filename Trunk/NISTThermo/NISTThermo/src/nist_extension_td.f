c
c*********************************************************************
c*********************************************************************
c***                                                               ***
c***            " NIST Extension Temperature and Density "         ***
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
c***                  subroutine TDSSOUND                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDSSOUND
c  
c  Flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDSSOUND (t,D,z,w,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDSSOUND
c     dll_export TDSSOUND
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
 1006   format ('[TDSSOUND warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDSSOUND error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDSSOUND error 220] ',a235,a1)
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
 1222     format ('[TDSSOUND error 222] dew point calculation ',
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
 1221       format ('[TDSSOUND error 221] bubble point calculation ',
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
c1225       format (1x,' TDSSOUND--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          w=xnotd     !Cp,w not defined for 2-phase states
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
      end                                         !subroutine TDSSOUND
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDDPDT2                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDDPDT2:
c  
c  Flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDDPDT2 (t,D,z,d2pdt2,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDDPDT2
c     dll_export TDDPDT2
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
 1006   format ('[TDDPDT2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDDPDT2 error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDDPDT2 error 220] ',a235,a1)
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
 1222     format ('[TDDPDT2 error 222] dew point calculation ',
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
 1221       format ('[TDDPDT2 error 221] bubble point calculation ',
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
c1225       format (1x,' TDDPDT2--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          d2pdt2=xnotd        !not defined for 2-phase states
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
      end                                          !subroutine TDDPDT2
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDDPDTD                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDDPDTD:
c  
c  Flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDDPDTD (t,D,z,d2pdtrho,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDDPDTD
c     dll_export TDDPDTD
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
 1006   format ('[TDDPDTD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDDPDTD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDDPDTD error 220] ',a235,a1)
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
 1222     format ('[TDDPDTD error 222] dew point calculation ',
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
 1221       format ('[TDDPDTD error 221] bubble point calculation ',
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
c1225       format (1x,' TDDPDTD--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          d2pdtrho=xnotd        !not defined for 2-phase states
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
      end                                          !subroutine TDDPDTD
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDDHDTCD                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDDHDTCD
c  
c  Flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDDHDTCD (t,D,z,dhdt_d,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDDHDTCD
c     dll_export TDDHDTCD
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
 1006   format ('[TDDHDTCD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDDHDTCD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDDHDTCD error 220] ',a235,a1)
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
 1222     format ('[TDDHDTCD error 222] dew point calculation ',
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
 1221       format ('[TDDHDTCD error 221] bubble point calculation ',
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
c1225       format (1x,' TDDHDTCD--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          dhdt_d = xnotd        !not defined for 2-phase states
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
      end                                         !subroutine TDDHDTCD
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDDHDTCP                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDDHDTCP
c  
c  Flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDDHDTCP (t,D,z,dhdt_p,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDDHDTCP
c     dll_export TDDHDTCP
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
 1006   format ('[TDDHDTCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDDHDTCP error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDDHDTCP error 220] ',a235,a1)
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
 1222     format ('[TDDHDTCP error 222] dew point calculation ',
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
 1221       format ('[TDDHDTCP error 221] bubble point calculation ',
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
c1225       format (1x,' TDDHDTCP--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          dhdt_p = xnotd        !not defined for 2-phase states
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
      end                                         !subroutine TDDHDTCP
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDDHDDCT                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDDHDDCT
c  
c  Flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDDHDDCT (t,D,z,dhdd_t,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDDHDDCT
c     dll_export TDDHDDCT
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
 1006   format ('[TDDHDDCT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDDHDDCT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDDHDDCT error 220] ',a235,a1)
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
 1222     format ('[TDDHDDCT error 222] dew point calculation ',
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
 1221       format ('[TDDHDDCT error 221] bubble point calculation ',
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
c1225       format (1x,' TDDHDDCT--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          dhdd_t = xnotd        !not defined for 2-phase states
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
      end                                         !subroutine TDDHDDCT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDDHDDCP                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDDHDDCP
c  
c  Flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDDHDDCP (t,D,z,dhdd_p,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDDHDDCP
c     dll_export TDDHDDCP
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
 1006   format ('[TDDHDDCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDDHDDCP error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDDHDDCP error 220] ',a235,a1)
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
 1222     format ('[TDDHDDCP error 222] dew point calculation ',
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
 1221       format ('[TDDHDDCP error 221] bubble point calculation ',
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
c1225       format (1x,' TDDHDDCP--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          dhdd_p = xnotd        !not defined for 2-phase states
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
      end                                         !subroutine TDDHDDCP
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDDHDPCT                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDDHDPCT
c  
c  Flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDDHDPCT (t,D,z,dhdp_t,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDDHDPCT
c     dll_export TDDHDPCT
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
 1006   format ('[TDDHDPCT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDDHDPCT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDDHDPCT error 220] ',a235,a1)
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
 1222     format ('[TDDHDPCT error 222] dew point calculation ',
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
 1221       format ('[TDDHDPCT error 221] bubble point calculation ',
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
c1225       format (1x,' TDDHDPCT--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          dhdp_t = xnotd        !not defined for 2-phase states
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
      end                                         !subroutine TDDHDPCT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine TDDHDPCD                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine TDDHDPCD
c  
c  Flash calculation given temperature, bulk density,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine TDDHDPCD (t,D,z,dhdp_d,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
c
cDEC$ ATTRIBUTES DLLEXPORT :: TDDHDPCD
c     dll_export TDDHDPCD
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
 1006   format ('[TDDHDPCD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[TDDHDPCD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=220
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[TDDHDPCD error 220] ',a235,a1)
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
 1222     format ('[TDDHDPCD error 222] dew point calculation ',
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
 1221       format ('[TDDHDPCD error 221] bubble point calculation ',
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
c1225       format (1x,' TDDHDPCD--TDFL2 return   p,q,x,y = ',8f12.7)
          end if
c  compute remaining properties for 2-phase states
          call PRESS (t,Dl,x,pl)
          call PRESS (t,Dv,y,pv)
          if (nc.ge.2) p=0.5d0*(pl+pv)
          dhdp_d = xnotd        !not defined for 2-phase states
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
      end                                         !subroutine TDDHDPCD
