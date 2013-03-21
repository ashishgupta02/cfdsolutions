c
c*********************************************************************
c*********************************************************************
c***                                                               ***
c***            " NIST Extension Pressure and Density "            ***
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
c***                  subroutine PDSSOUND                          ***
c***                                                               ***
c*********************************************************************
c
c           Copyright 2013 CSEWorks INC All Rights Reserved
c                      PROPRIETARY INFORMATION
c                   For CSEWorks INC Use Only.
c
c
c  Subroutine PDSSOUND
c  
c  Flash calculation given density, pressure,
c  and bulk composition
c
c---------------------------------------------------------------------
      subroutine PDSSOUND (p,D,z,w,ierr,herr)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      implicit logical (l)
c
cDEC$ ATTRIBUTES DLLEXPORT :: PDSSOUND
c     dll_export PDSSOUND
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
 1008   format ('[PDSSOUND error 210] ',a235,a1)
        call ERRMSG (ierr,herr)
        RETURN
      end if
c
      t=-1
      call LIMITX ('EOS',t,D,p,z,tmin,tmax,rhomax,pmax,ierrl,herrl)
      if (ierrl.lt.0) then
c  one or inputs are outside limits--if just a warning proceed w/ calc
        write (herr,1006) ierrl,herrl(1:234),hnull
 1006   format ('[PDSSOUND warning',i3,'] ',a234,a1)
        call ERRMSG (ierrl,herr)
c  if error (as opposed to warning) set output density to zero and return
      else if (ierrl.gt.0) then
        write (herr,1007) ierrl,herrl(1:236),hnull
 1007   format ('[PDSSOUND error',i3,'] ',a236,a1)
        call ERRMSG (ierrl,herr)
        ierr=ierrl
        RETURN
      end if
c
      if (p.ge.pc) then
c  supercritical state (x = y = z as set above)
c       write (*,*) ' PDSSOUND--supercritical'
        call PDFL1 (p,D,z,t,ierr,herr2)
        if (ierr.ne.0) then
          ierr=215
          write (herr,1215) herr2(1:183),hnull
 1215   format('[PDSSOUND error 215] supercritical density iteration ',
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
 1212     format ('[PDSSOUND error 212] dew point calculation ',
     &            'did not converge:  ',a193,a1)
          call ERRMSG (ierr,herr)
          RETURN
        elseif (ierr.ne.0) then
c  do not exit yet, check if valid liquid state
          iflag=1
        end if
        if (D.le.Dvdew+1.0d-12) then
c  single-phase vapor (x = y = z as set above)
c         write (*,*) ' PDSSOUND--single-phase vapor'
          call PDFL1 (p,D,z,t,ierr,herr2)
          if (ierr.ne.0) then
            ierr=214
            write (herr,1214) herr2(1:190),hnull
 1214       format ('[PDSSOUND error 214] vapor density iteration ',
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
 1211       format ('[PDSSOUND error 211] bubble point calculation ',
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
 1213         format ('[PDSSOUND error 213] liquid density iteration ',
     &                'did not converge:  ',a189,a1)
              call ERRMSG (ierr,herr)
              t=0.0d0
              RETURN
            end if
            call THERM (t,D,z,ptherm,e,h,s,cv,cp,w,hjt)
          else
c
c  two-phase pure fluid
            w=xnotd     !Cp,w not defined for 2-phase states
            RETURN
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
 2006   format ('[PDSSOUND warning',i3,'] ',a234,a1)
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
      end                                         !subroutine PDSSOUND
c
c
