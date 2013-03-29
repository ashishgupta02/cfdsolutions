c
c*********************************************************************
c*********************************************************************
c***                                                               ***
c***      " NIST Extension Pressure and Density Multiphase "       ***
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
c***                  subroutine PDEPDFLSH                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEPDFLSH: Extended PDFLSH
c  
c  Extended flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEPDFLSH (p,D,z,t,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEPDFLSH
c     dll_export PDEPDFLSH
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      Dl=0.d0
      Dv=0.d0
      e=0.d0
      h=0.d0
      s=0.d0
      cv=0.d0
      cp=0.d0
      w=0.d0
      q=999.0d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEPDFLSH error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEPDFLSH warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEPDFLSH error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEPDFLSH--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEPDFLSH error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        Dl=D
        Dv=D
        call THERM (t,D,z,ptherm,e,h,s,cv,cp,w,hjt)
      else
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEPDFLSH error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEPDFLSH--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEPDFLSH error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          Dl=D
          Dv=D
          call THERM (t,D,z,ptherm,e,h,s,cv,cp,w,hjt)
          q=998.d0
          if (tdew.gt.0 .and. ianc(icomp).eq.0) then
            call ENTHAL (tdew,Dldew,xdew,hldew)
            call ENTHAL (tdew,Dvdew,z,hvdew)
            q=(h-hldew)/(hvdew-hldew)
          endif
        else
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEPDFLSH error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEPDFLSH error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            Dl=D
            Dv=D
            call THERM (t,D,z,ptherm,e,h,s,cv,cp,w,hjt)
c  compute quality based on volumes
            q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
c           write (*,*) ' PDEPDFLSH--comp liq q by volumes:  ',q
          else
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEPDFLSH--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call THERM (t,Dl,x,ptherm,el,hl,sl,cvl,cpl,wl,hjt)
            call THERM (t,Dv,y,ptherm,ev,hv,sv,cvv,cpv,wv,hjt)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            e=alpha*el+q*ev
            h=alpha*hl+q*hv
            s=alpha*sl+q*sv
c  compute the harmonical average for 2-phase states
            tmp=(q*Dv+alpha*Dl)
            w=tmp*(q/(Dv*wv*wv) + alpha/(Dl*wl*wl))
            w=SQRT(1.0d0/w);
            cp=tmp*(q/(Dv*cpv*cpv) + alpha/(Dl*cpl*cpl))
            cp=SQRT(1.0d0/cp);
            cv=tmp*(q/(Dv*cvv*cvv) + alpha/(Dl*cvl*cvl))
            cv=SQRT(1.0d0/cv);
          end if
        end if
      end if
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEPDFLSH warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDEPDFLSH
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEPDFLSH2                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEPDFLSH2: Extended PDFLSH with derivatives
c  
c  Extended flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEPDFLSH2 (p,D,z,x,y,thermv,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEPDFLSH2
c     dll_export PDEPDFLSH2
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension thermv(29)
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      do i=1,29
        thermv(i) = 0.d0
      enddo
      ierr = 0
      herr = ' '
      t    = 0.d0
      Dl   = 0.d0
      Dv   = 0.d0
      q    = 999.0d0
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
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        thermv(1) = t
        thermv(4) = q
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEPDFLSH2 error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        thermv(1) = t
        thermv(4) = q
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEPDFLSH2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEPDFLSH2 error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        thermv(1) = t
        thermv(4) = q
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEPDFLSH2--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215  format('[PDEPDFLSH2 error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          thermv(1) = t
          thermv(4) = q
          RETURN
        end if
        Dl=D
        Dv=D
        call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEPDFLSH2 error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          thermv(1) = t
          thermv(4) = q
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEPDFLSH2--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEPDFLSH2 error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            thermv(1) = t
            thermv(4) = q
            RETURN
          end if
          Dl=D
          Dv=D
          call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,
     &      beta,dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          q=998.d0
          if (tdew.gt.0 .and. ianc(icomp).eq.0) then
            call ENTHAL (tdew,Dldew,xdew,hldew)
            call ENTHAL (tdew,Dvdew,z,hvdew)
            q=(h-hldew)/(hvdew-hldew)
          endif
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEPDFLSH2 error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            thermv(1) = t
            thermv(4) = q
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213        format ('[PDEPDFLSH2 error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              thermv(1) = t
              thermv(4) = q
              RETURN
            end if
            Dl=D
            Dv=D
            call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,
     &        beta,dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &        spare3,spare4)
            call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
c  compute quality based on volumes
            q=(1.0d0/D-1.0d0/Dlbub)/(1.0d0/Dvdew-1.0d0/Dlbub)
c           write (*,*) ' PDEPDFLSH2--comp liq q by volumes:  ',q
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEPDFLSH2--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                thermv(1) = t
                thermv(4) = q
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                thermv(1) = t
                thermv(2) = Dl
                thermv(3) = Dv
                thermv(4) = q
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call THERM2 (t,Dl,x,ptherm,el,hl,sl,cvl,cpl,wl,zcl,hjtl,
     &        al,gl,xkappal,betal,dpdrhol,d2pdd2l,dptl,drhodtl,drhodpl,
     &        d2pdt2l,d2pdtrhol,spare3,spare4)
c
            call THERM2 (t,Dv,y,ptherm,ev,hv,sv,cvv,cpv,wv,zcv,hjtv,
     &        av,gv,xkappav,betav,dpdrhov,d2pdd2v,dptv,drhodtv,drhodpv,
     &        d2pdt2v,d2pdtrhov,spare3,spare4)
c
            call DHD1 (t,Dl,x,dhdt_dl,dhdt_pl,dhdd_tl,dhdd_pl,
     &      dhdp_tl,dhdp_dl)
c
            call DHD1 (t,Dv,y,dhdt_dv,dhdt_pv,dhdd_tv,dhdd_pv,
     &      dhdp_tv,dhdp_dv)
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
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  Finally return the thermodynamic properties
      thermv(1)  = t
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
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEPDFLSH2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                       !subroutine PDEPDFLSH2
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDETEMP                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDETEMP
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDETEMP (p,D,z,t,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDETEMP
c     dll_export PDETEMP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDETEMP error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDETEMP warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDETEMP error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDETEMP--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDETEMP error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
      else
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDETEMP error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDETEMP--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDETEMP error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
        else
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDETEMP error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDETEMP error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
          else
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
          end if
        end if
      end if
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDETEMP warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                          !subroutine PDETEMP
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDESSOUND                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDESSOUND
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDESSOUND (p,D,z,w,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDESSOUND
c     dll_export PDESSOUND
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      w=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDESSOUND error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDESSOUND warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDESSOUND error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDESSOUND--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDESSOUND error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        Dl=D
        Dv=D
        call THERM (t,D,z,ptherm,e,h,s,cv,cp,w,hjt)
      else
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDESSOUND error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDESSOUND--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDESSOUND error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call THERM (t,D,z,ptherm,e,h,s,cv,cp,w,hjt)
        else
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDESSOUND error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDESSOUND error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call THERM (t,D,z,ptherm,e,h,s,cv,cp,w,hjt)
          else
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDESSOUND--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call THERM (t,Dl,x,ptherm,el,hl,sl,cvl,cpl,wl,hjt)
            call THERM (t,Dv,y,ptherm,ev,hv,sv,cvv,cpv,wv,hjt)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
c  compute the harmonical average for 2-phase states
            tmp=(q*Dv+alpha*Dl)
            w=tmp*(q/(Dv*wv*wv) + alpha/(Dl*wl*wl))
            w=SQRT(1.0d0/w);
          end if
        end if
      end if
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDESSOUND warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDESSOUND
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEENTRO                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEENTRO
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEENTRO (p,D,z,s,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEENTRO
c     dll_export PDEENTRO
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      s=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEENTRO error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEENTRO warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEENTRO error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEENTRO--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEENTRO error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call ENTRO (t,D,z,s)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEENTRO error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEENTRO--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEENTRO error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call ENTRO (t,D,z,s)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEENTRO error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEENTRO error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call ENTRO (t,D,z,s)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEENTRO--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call ENTRO (t,Dl,x,sl)
            call ENTRO (t,Dv,y,sv)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            s=alpha*sl+q*sv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEENTRO warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                         !subroutine PDEENTRO
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEENTHAL                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEENTHAL
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEENTHAL (p,D,z,h,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEENTHAL
c     dll_export PDEENTHAL
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      h=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEENTHAL error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEENTHAL warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEENTHAL error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEENTHAL--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEENTHAL error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call ENTHAL (t,D,z,h)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEENTHAL error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEENTHAL--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEENTHAL error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call ENTHAL (t,D,z,h)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEENTHAL error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEENTHAL error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call ENTHAL (t,D,z,h)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEENTHAL--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call ENTHAL (t,Dl,x,hl)
            call ENTHAL (t,Dv,y,hv)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            h=alpha*hl+q*hv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEENTHAL warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDEENTHAL
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEENERGY                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEENERGY
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEENERGY (p,D,z,e,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEENERGY
c     dll_export PDEENERGY
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      e=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEENERGY error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEENERGY warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEENERGY error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEENERGY--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEENERGY error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call ENERGY (t,D,z,e)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEENERGY error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEENERGY--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEENERGY error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call ENERGY (t,D,z,e)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEENERGY error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEENERGY error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call ENERGY (t,D,z,e)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEENERGY--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call ENERGY (t,Dl,x,el)
            call ENERGY (t,Dv,y,ev)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            e=alpha*el+q*ev
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEENERGY warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDEENERGY
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDECVCP                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDECVCP
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDECVCP (p,D,z,cv,cp,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDECVCP
c     dll_export PDECVCP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      cv=0.d0
      cp=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDECVCP error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDECVCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDECVCP error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDECVCP--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDECVCP error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call CVCP (t,D,z,cv,cp)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDECVCP error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDECVCP--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDECVCP error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call CVCP (t,D,z,cv,cp)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDECVCP error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDECVCP error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call CVCP (t,D,z,cv,cp)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDECVCP--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call CVCP (t,Dl,x,cvl,cpl)
            call CVCP (t,Dv,y,cvv,cpv)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            tmp=(q*Dv+alpha*Dl)
            cp=tmp*(q/(Dv*cpv*cpv) + alpha/(Dl*cpl*cpl))
            cp=SQRT(1.0d0/cp);
            cv=tmp*(q/(Dv*cvv*cvv) + alpha/(Dl*cvl*cvl))
            cv=SQRT(1.0d0/cv);
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDECVCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                          !subroutine PDECVCP
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDPDD                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDPDD
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDPDD (p,D,z,dpdrho,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDPDD
c     dll_export PDEDPDD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      dpdrho=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDPDD error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDPDD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDPDD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDPDD--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDPDD error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DPDD (t,D,z,dpdrho)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDPDD error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDPDD--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDPDD error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DPDD (t,D,z,dpdrho)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDPDD error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDPDD error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DPDD (t,D,z,dpdrho)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDPDD--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DPDD (t,Dl,x,dpdrhol)
            call DPDD (t,Dv,y,dpdrhov)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            dpdrho=alpha*dpdrhol+q*dpdrhov
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDPDD warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                          !subroutine PDEDPDD
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDPDT                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDPDT
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDPDT (p,D,z,dpt,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDPDT
c     dll_export PDEDPDT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      dpt=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDPDT error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDPDT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDPDT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDPDT--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDPDT error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DPDT (t,D,z,dpt)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDPDT error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDPDT--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDPDT error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DPDT (t,D,z,dpt)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDPDT error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDPDT error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DPDT (t,D,z,dpt)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDPDT--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DPDT (t,Dl,x,dptl)
            call DPDT (t,Dv,y,dptv)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            dpt=alpha*dptl+q*dptv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDPDT warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                          !subroutine PDEDPDT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDDDT                           ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDDDT
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDDDT (p,D,z,drhodt,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDDDT
c     dll_export PDEDDDT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      drhodt=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDDDT error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDDDT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDDDT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDDDT--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDDDT error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DDDT (t,D,z,drhodt)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDDDT error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDDDT--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDDDT error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DDDT (t,D,z,drhodt)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDDDT error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDDDT error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DDDT (t,D,z,drhodt)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDDDT--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DDDT (t,Dl,x,drhodtl)
            call DDDT (t,Dv,y,drhodtv)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            drhodt=alpha*drhodtl+q*drhodtv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDDDT warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                          !subroutine PDEDDDT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDPDD2                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDPDD2
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDPDD2 (p,D,z,d2pdd2,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDPDD2
c     dll_export PDEDPDD2
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      d2pdd2=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDPDD2 error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDPDD2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDPDD2 error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDPDD2--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDPDD2 error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DPDD2 (t,D,z,d2pdd2)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDPDD2 error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDPDD2--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDPDD2 error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DPDD2 (t,D,z,d2pdd2)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDPDD2 error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDPDD2 error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DPDD2 (t,D,z,d2pdd2)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDPDD2--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DPDD2 (t,Dl,x,d2pdd2l)
            call DPDD2 (t,Dv,y,d2pdd2v)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            d2pdd2=alpha*d2pdd2l+q*d2pdd2v
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDPDD2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                         !subroutine PDEDPDD2
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDPDT2                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDPDT2
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDPDT2 (p,D,z,d2pdt2,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDPDT2
c     dll_export PDEDPDT2
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      d2pdt2=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDPDT2 error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDPDT2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDPDT2 error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDPDT2--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDPDT2 error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDPDT2 error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDPDT2--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDPDT2 error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,
     &      beta,dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDPDT2 error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDPDT2 error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,
     &      beta,dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDPDT2--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call THERM2 (t,Dl,x,pl,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2l,d2pdtrho,
     &      spare3,spare4)
            call THERM2 (t,Dv,y,pv,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2v,d2pdtrho,
     &      spare3,spare4)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            d2pdt2=alpha*d2pdt2l+q*d2pdt2v
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDPDT2 warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                         !subroutine PDEDPDT2
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDPDTD                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDPDTD
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDPDTD (p,D,z,d2pdtrho,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDPDTD
c     dll_export PDEDPDTD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      d2pdtrho=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDPDTD error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDPDTD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDPDTD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDPDTD--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDPDTD error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDPDTD error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDPDTD--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDPDTD error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,
     &      beta,dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDPDTD error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDPDTD error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call THERM2 (t,D,z,ptherm,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,
     &      beta,dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrho,
     &      spare3,spare4)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDPDTD--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call THERM2 (t,Dl,x,pl,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrhol,
     &      spare3,spare4)
            call THERM2 (t,Dv,y,pv,e,h,s,cv,cp,w,zc,hjt,a,g,xkappa,beta,
     &      dpdrho,d2pdd2,dpt,drhodt,drhodp,d2pdt2,d2pdtrhov,
     &      spare3,spare4)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            d2pdtrho=alpha*d2pdtrhol+q*d2pdtrhov
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDPDTD warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                         !subroutine PDEDPDTD
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDHDTCD                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDHDTCD
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDHDTCD (p,D,z,dhdt_d,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDHDTCD
c     dll_export PDEDHDTCD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      dhdt_d=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDHDTCD error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDHDTCD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDHDTCD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDHDTCD--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDHDTCD error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDHDTCD error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDHDTCD--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDHDTCD error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDHDTCD error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDHDTCD error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDHDTCD--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DHD1 (t,Dl,x,dhdt_dl,dhdt_p,dhdd_t,dhdd_p,
     &      dhdp_t,dhdp_d)
            call DHD1 (t,Dv,y,dhdt_dv,dhdt_p,dhdd_t,dhdd_p,
     &      dhdp_t,dhdp_d)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            dhdt_d=alpha*dhdt_dl+q*dhdt_dv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDHDTCD warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDEDHDTCD
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDHDTCP                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDHDTCP
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDHDTCP (p,D,z,dhdt_p,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDHDTCP
c     dll_export PDEDHDTCP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      dhdt_p=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDHDTCP error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDHDTCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDHDTCP error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDHDTCP--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDHDTCP error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDHDTCP error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDHDTCP--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDHDTCP error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDHDTCP error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDHDTCP error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDHDTCP--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DHD1 (t,Dl,x,dhdt_d,dhdt_pl,dhdd_t,dhdd_p,
     &      dhdp_t,dhdp_d)
            call DHD1 (t,Dv,y,dhdt_d,dhdt_pv,dhdd_t,dhdd_p,
     &      dhdp_t,dhdp_d)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            dhdt_p=alpha*dhdt_pl+q*dhdt_pv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDHDTCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDEDHDTCP
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDHDDCT                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDHDDCT
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDHDDCT (p,D,z,dhdd_t,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDHDDCT
c     dll_export PDEDHDDCT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      dhdd_t=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDHDDCT error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDHDDCT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDHDDCT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDHDDCT--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDHDDCT error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDHDDCT error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDHDDCT--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDHDDCT error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDHDDCT error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDHDDCT error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDHDDCT--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DHD1 (t,Dl,x,dhdt_d,dhdt_p,dhdd_tl,dhdd_p,
     &      dhdp_t,dhdp_d)
            call DHD1 (t,Dv,y,dhdt_d,dhdt_p,dhdd_tv,dhdd_p,
     &      dhdp_t,dhdp_d)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            dhdd_t=alpha*dhdd_tl+q*dhdd_tv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDHDDCT warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDEDHDDCT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDHDDCP                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDHDDCP
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDHDDCP (p,D,z,dhdd_p,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDHDDCP
c     dll_export PDEDHDDCP
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      dhdd_p=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDHDDCP error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDHDDCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDHDDCP error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDHDDCP--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDHDDCP error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDHDDCP error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDHDDCP--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDHDDCP error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDHDDCP error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDHDDCP error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDHDDCP--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DHD1 (t,Dl,x,dhdt_d,dhdt_p,dhdd_t,dhdd_pl,
     &      dhdp_t,dhdp_d)
            call DHD1 (t,Dv,y,dhdt_d,dhdt_p,dhdd_t,dhdd_pv,
     &      dhdp_t,dhdp_d)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            dhdd_p=alpha*dhdd_pl+q*dhdd_pv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDHDDCP warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDEDHDDCP
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDHDPCT                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDHDPCT
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDHDPCT (p,D,z,dhdp_t,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDHDPCT
c     dll_export PDEDHDPCT
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      dhdp_t=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDHDPCT error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDHDPCT warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDHDPCT error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDHDPCT--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDHDPCT error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDHDPCT error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDHDPCT--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDHDPCT error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDHDPCT error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDHDPCT error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDHDPCT--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DHD1 (t,Dl,x,dhdt_d,dhdt_p,dhdd_t,dhdd_p,
     &      dhdp_tl,dhdp_d)
            call DHD1 (t,Dv,y,dhdt_d,dhdt_p,dhdd_t,dhdd_p,
     &      dhdp_tv,dhdp_d)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            dhdp_t=alpha*dhdp_tl+q*dhdp_tv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDHDPCT warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDEDHDPCT
c
c
c*********************************************************************
c***                                                               ***
c***                  subroutine PDEDHDPCD                         ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDEDHDPCD
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDEDHDPCD (p,D,z,dhdp_d,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDEDHDPCD
c     dll_export PDEDHDPCD
c
      parameter (ncmax=20)        !max number of components in mixture
      parameter (nrefmx=10)       !max number of fluids for transport E
      parameter (n0=-ncmax-nrefmx,nx=ncmax)
      character*1 htab,hnull
      character*255 herr,herr2,herrl
      dimension z(ncmax),x(ncmax),y(ncmax)
      dimension xdew(ncmax),ydew(ncmax),ybub(ncmax),xbub(ncmax)
      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
      common /NCOMP/ nc,ic
      common /HCHAR/ htab,hnull
c  flags to GUI indicating 'not applicable', '2-phase', etc.
      common /FLAGS/ xnota,x2ph,xsubc,xsuph,xsupc,xinf,x7,xnotd,xnotc
      common /FLAGS2/ iamwat,ianc(0:ncmax),iwat
      common /CCON/ tcrit(n0:nx),pcrit(n0:nx),Dcrit(n0:nx),Zcrit(n0:nx),
     &              ttp(n0:nx),ptp(n0:nx),dtp(n0:nx),dtpv(n0:nx),
     &              tnbp(n0:nx),dnbpl(n0:nx),dnbpv(n0:nx),
     &              wm(n0:nx),accen(n0:nx),dipole(n0:nx),Reos(n0:nx)
      common /prnterr/ iprnterr
c
      call ISPURE (z,icomp)
      ierr=0
      herr=' '
      t=0.d0
      dhdp_d=0.d0
c
c  initialize output liquid and vapor compositions to input values
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

      if (p.lt.1.0d-14 .or. d.lt.1.0d-14) then
c  ideal-gas
        if (D.gt.0.d0) t=p/(R*D)
        RETURN
      end if
c
c  check that input conditions are within limits
c
      call CRITP (z,tc,pc,rhoc,ierr,herr2)
      if (ierr.gt.0) then
        ierr=210
        write (herr,1008) herr2(1:235),hnull
 1008   format ('[PDEDHDPCD error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDEDHDPCD warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDEDHDPCD error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then                                         !if 0
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDEDHDPCD--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDEDHDPCD error 215] supercritical density iteration ',
     &            'did not converge:  ',a183,a1)
          call ERRMSG (ierr,herr)
          RETURN
        end if
        call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
      else                                                     !else 0
c  subcritical state--call saturation routine to determine liq or vap
        i=iprnterr
        iprnterr=0
        ierr=0
        Dvdew=D
        Dldew=D
        tdew=0
c  do not call SATP if P<triple point pressure for pure fluids
        if (icomp.eq.0 .or. p.gt.ptp(icomp) .or. ianc(icomp).ne.0) then
          call SATP (p,z,2,tdew,Dldew,Dvdew,xdew,ydew,ierr,herr2)
        endif
        iprnterr=i
        iflag=0
        if (ierr.ne.0 .and. icomp.ne.0) then
          ierr=212
          write (herr,1212) herr2(1:193),hnull
 1212     format ('[PDEDHDPCD error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then                            !if 1
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDEDHDPCD--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDEDHDPCD error 214] vapor density iteration ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            t=0.0d0
            RETURN
          end if
          call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
        else                                                  !else 1
c  mixture:  calculate bubble point to determine if liquid or 2-phase
          if (icomp.ne.0 .and. ianc(icomp).eq.0) then
            Dlbub=Dldew
            Dvbub=Dvdew
            tbub=tdew
          else
            call SATP (p,z,1,tbub,Dlbub,Dvbub,xbub,ybub,ierr,herr2)
          endif
          if (ierr.ne.0) then
            ierr=211
            write (herr,1211) herr2(1:190),hnull
 1211       format ('[PDEDHDPCD error 211] bubble point calculation ',
     &              'did not converge:  ',a190,a1)
            call ERRMSG (ierr,herr)
            RETURN
          end if
          if (D.gt.Dlbub) then                                  !if 2
c  mixture single-phase liquid
c           write (*,*) ' TPFLSH--mixture single-phase liquid'
            call PDFL1 (p,D,z,t,ierr,herr2)
            if (ierr.ne.0) then
              ierr=213
              write (herr,1213) herr2(1:189),hnull
 1213         format ('[PDEDHDPCD error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call DHD1 (t,D,z,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
          else                                                !else 2
c
c  two-phase pure fluid
            if (icomp.ne.0) then
              Dl=Dlbub
              Dv=Dvdew
              q=(1.0d0/D-1.0d0/Dl)/(1.0d0/Dv-1.0d0/Dl)
              t=(1.d0-q)*tbub+q*tdew
              if (ianc(icomp).eq.1) call PTANC (t,p,q,d,'D',Dl,Dv)
c  two-phase mixture
c             write (*,*) ' PDEDHDPCD--two-phase mixture'
c
c  generate initial guesses for x and y by interpolating input pressure
c  with dew and bubble point pressures
c
            else
              if (iflag.eq.1) then
c  dew point did not converge, location of state point is either two-phase
c  or vapor phase
                ierr=212
                write (herr,1212) herr2(1:193),hnull
                call ERRMSG (ierr,herr)
                t=0.0d0
                RETURN
              endif
              ksat=1   !bubble and dew point data provided to PDFL2
              call PDFL2 (p,D,z,ksat,tbub,tdew,Dlbub,Dvdew,ybub,xdew,
     &                    t,Dl,Dv,x,y,q,ierr,herr)
              if (ierr.ne.0) then
c  two-phase iteration did not converge--error message written by PDFL2
                q=999.d0     !quality undefined
                RETURN
              end if
c             write (*,1225) p,q,(x(i),i=1,3),(y(i),i=1,3)
c1225         format (1x,' TDFLSH--TDFL2 return   p,q,x,y = ',8f12.7)
            endif
c  compute remaining properties for 2-phase states
            call DHD1 (t,Dl,x,dhdt_d,dhdt_p,dhdd_t,dhdd_p,
     &      dhdp_t,dhdp_dl)
            call DHD1 (t,Dv,y,dhdt_d,dhdt_p,dhdd_t,dhdd_p,
     &      dhdp_t,dhdp_dv)
c  bulk properties are weighted average of liquid and vapor phases
            alpha=1.0d0-q               !alpha is liq fraction,
            dhdp_d=alpha*dhdp_dl+q*dhdp_dv
          end if                                            !end if 2
        end if                                              !end if 1
      end if                                                !end if 0
c
c  if limits check resulted in warning (as opposed to error) return
c  that message; do this again in case intermediate iteration did not
c  converge (thereby overwriting any warning message from LIMITX)
c
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        ierr=ierrl
        write (herr,2006) ierrl,herrl(1:234),hnull
 2006   format ('[PDEDHDPCD warning',i3,'] ',a234,a1)
        call ERRMSG (ierr,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
        ierr=ierrl
        call ERRMSG (ierrl,herr)
        RETURN
      end if
c
      RETURN
      end                                        !subroutine PDEDHDPCD
c
c
