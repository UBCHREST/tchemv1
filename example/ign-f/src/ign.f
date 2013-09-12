      program ign

      implicit none 
      include 'global.par'
      real*8 runical
      integer liwmax, lrwmax
      parameter (runical = 1.9872065009560d0)
      parameter (liwmax = 30+neqmax, lrwmax = 22+9*neqmax+2*neqmax**2)

      integer  i, iter, itermax, ispec
      integer  iopt, iout, istate, itask, itol, iwork(liwmax),liw,
     &         lrw, mf, neq, ns, lmech, lthrm, mktab, outfreq
      real*8   atol(neqmax), rtol, stol, rwork(lrwmax), 
     &         t, tout, y(neqmax), ymfr(nsmax)
      real*8   deltaTemp, deltatsta, deltatmax, deltatmin, deltat, tmax
      real*8   Tempini, Told, Xo2, Xch4, Xar, Xn2, dtemptab
      real*8   Temp_id, Temp_m2, Temp_m1, time_m2, time_m1, time_id1, 
     &         der2num, time_id2
      integer  foundid2
      integer  reacid,posid
      real*8   preexp,tempexp,acten

      character mechfile*100, thermofile*100
      character*18 specnames(nsmax)
      integer ipar(1)
      real*8  rpar(1)

      integer      nspecin
      real*8       smsfrin(neqmax)
      character*18 snamein(neqmax)

      external fex, jex
      external tcgetml2ms, tcinitchem, tcgetms2hmixms, tcgetms2cc,
     &         tcgetsnames, tcgetarhenfor, tcchgarhenfor, tcreset,
     &         tcgetnspec
      integer tcgetml2ms, tcinitchem, tcgetms2hmixms, tcgetms2cc, 
     &        tcgetsnames, tcgetarhenfor, tcchgarhenfor, tcgetnspec, 
     &        ierr

c-Max temp change, time steps (start, max, min), and end time
      deltaTemp = 1.0d0
      deltatsta = 1.0d-5  
      deltatmax = 1.0d-3
      deltatmin = 1.0d-11
      tmax      = 1.0d1
      itermax   = 1000000
      outfreq   = 10
      Tempini   = 1000.d0
      Temp_id   = 1500.d0

      mechfile   = 'chem.inp'
      thermofile = 'therm.dat'
      mktab      = 0
      dtemptab   = 0.1d0

c-read from input file
      call initrun(itermax, outfreq, mktab,deltatsta, deltatmax, 
     &             tmax, Tempini, Temp_id, deltaTemp, smsfrin, 
     &             snamein, nspecin, mechfile, thermofile)
      write(*,*)'No. of species in the setup file    : ',nspecin

c-initialize tchem
      lmech = len_trim(mechfile  )
      lthrm = len_trim(thermofile)
      ierr = tcinitchem(mechfile,lmech,thermofile,lthrm, 
     &                  mktab, dtemptab)

c-initialize fresh mixture
      ns = tcgetnspec() 
      write(*,*)'No. of species in the kinetic model : ',ns
      neq = ns + 1
      do i = 1,neq
         y(i) = 0.d0
      end do

      ierr = tcgetsnames( ns, specnames )
      do ispec = 1,nspecin
         do i = 1,ns
            if (specnames(i).eq.snamein(ispec)) then
               y(i+1) = smsfrin(ispec)
               write(*,'(A8,A,A4,I4)')' found ',specnames(i),' at ',i
               exit
            endif
         enddo
      enddo
      y(1) = Tempini 

c-transform mole fractions to mass fractions
      ierr = tcgetml2ms( y(2), ns, ymfr )
      do i=1,ns
        y(i+1)=ymfr(i)
      end do

#ifdef MODKINMECH
c-Print Arrhenius parameters for reaction 0
      reacid = 0
      posid  = 0
      ierr = tcgetarhenfor(reacid,posid,preexp)
      posid  = 1
      ierr = tcgetarhenfor(reacid,posid,tempexp)
      posid  = 2
      ierr = tcgetarhenfor(reacid,posid,acten)
      write(*,'(A,I4,1x,E20.12,1x,F9.4,1x,E20.12)')
     &      'Arrhenius parameters :',reacid,preexp,tempexp,acten

      reacid = 0
      posid  = 2
      acten = acten*1.4d0
      ierr = tcchgarhenfor(reacid,posid,acten) 
      posid  = 2
      ierr = tcgetarhenfor(reacid,posid,acten)
      write(*,'(A,I4,1x,E20.12,1x,F9.4,1x,E20.12)')
     &      'Arrhenius parameters :',reacid,preexp,tempexp,acten
#endif

c-Output header files and open solution files
#ifndef NO_OUTPUT
      open( unit=10, file='ignsol.hdr' )
      open( unit=11, file='ys.hdr'     )
      open( unit=12, file='cs.hdr'     )
      open( unit=13, file='h.hdr'      )
      write(10,'(A36)')'# 1: iteration, #2: time, #3: deltat'
      write(10,10)(i+3,specnames(i)(1:len_trim(specnames(i))),i=1,ns)
      write(11,'(A9)')'# 1: time'
      write(11,10)(i+1,specnames(i),i=1,ns)
      write(12,'(A9)')'# 1: time'
      write(12,12)(i+1,specnames(i),i=1,ns)
      write(13,'(A39)')'# 1: time, #2: specific enthalpy [J/kg]'
 10   format('#',i3,' ',a18,' Mass Fraction '      )
 12   format('#',i3,' ',a18,' Molar Concentration [kmol/m3]')
      close(10)
      close(11)
      close(12)
      close(13)

      open( unit=10, file='ignsol.dat' )
      open( unit=11, file='ys.out'     )
      open( unit=12, file='cs.out'     )
      open( unit=13, file='h.out'      )
#endif

c-dvode setup
      itol = 2
      rtol = 1.d-12
      stol = 1.d-20
      do i = 1,neq
        atol = max( rtol*y(i), stol )
      end do

      itask  = 1
      istate = 1
      iopt   = 0
      rwork(5) = deltatsta 
      rwork(6) = deltatmax
      rwork(7) = deltatmin
      iwork(5) = 5
      iwork(6) = 100000
      iwork(7) = 10
      lrw      = lrwmax
      liw      = liwmax
      mf       = 21  ! stiff, BDF method, user supplied Jacobian


c-Initialize time and temperature history
      foundid2 = 0
      Temp_m2 = Tempini
      Temp_m1 = Tempini
      time_m2 = 0.d0
      time_m1 = 0.d0

c-Start time advancement
      t    = 0.d0
      deltat = deltatsta
      iter = 0
      do while ( (t+deltat.le.tmax).and.(iter.le.itermax) )
         iter = iter+1
         tout = t+deltat
         Told = y(1)
         call dvode(fex, neq, y, t, tout, itol, rtol, atol, itask,
     &              istate, iopt, rwork, lrw, iwork, liw, jex, mf, 
     &              rpar, ipar)

c-calculate new time step:
c       - reduce/increase time step if \Delta T >< \Delta T_{set}
c       -  
         deltat = min(deltaTemp/max(abs(y(1)-Told),stol)*deltat,
     &                2.0*deltat
     &               )
         deltat = min(deltat,deltatmax)
         if ( tout+deltat.gt.tmax ) then
            deltat = max(tmax-tout,deltatmin)
         end if 
         t = tout

#ifndef NO_OUTPUT
c-Output solution
         if ( mod(iter,outfreq).eq.0 ) then
            write(6,20)  iter, tout, deltat, y(1)
20          format(' Iter = ',i7,'  t[s] = ',d14.6,'  dt[s] = ',d14.6,
     &            '   T[K] = ',d14.6)
            write(10,'(i10,1x,d20.12,1x,d20.12,1x,200(d20.12,1x))')
     &            iter, t, deltat, (y(i),i=1,neq)
            write(11,'(d20.12,1x,200(d20.12,1x))')
     &            t, (y(i),i=1,neq)
            ierr = tcgetms2cc( y, neq, ymfr ) ;
            write(12,'(d20.12,1x,200(d20.12,1x))')
     &            t, (ymfr(i),i=1,ns)
           ierr = tcgetms2hmixms( y, neq, ymfr(1) ) ;
           write(13,'(d20.12,1x,d20.12)') t, ymfr(1)

         end if ! done output solution
#endif

c-Ignition delay time - cross temperature threshold
         if ( (y(1)-Temp_id)*(Temp_m1-Temp_id).le.0.d0 ) then
            time_id1 = time_m1+(t-time_m1)
     &                *(Temp_id-Temp_m1)/(y(1)-Temp_m1)
         endif
c-Ignition delay time - second derivative
         if ( (foundid2.eq.0).and.(Temp_m1.gt.1300.0) ) then
            der2num = (y(1)   -Temp_m1)/(t      -time_m1)
     &               -(Temp_m1-Temp_m2)/(time_m1-time_m2)
            if ( der2num.le.0.d0 ) then
               foundid2 = 1
               time_id2 = time_m1
            endif
         endif

c-Cycle temperatures and times
         Temp_m2 = Temp_m1
         Temp_m1 = y(1)
         time_m2 = time_m1
         time_m1 = t

c-Need to recompute absolute tolerances
         do i = 1,neq
           atol = max( rtol*y(i), stol )
         end do

      end do ! end time advancement loop


#ifndef NO_OUTPUT
c-Output last solution
      write(6,20)  iter, tout, deltat, y(1)
      write(10,'(i10,1x,d20.12,1x,d20.12,1x,200(d20.12,1x))')
     &           iter, t, deltat, (y(i),i=1,neq)
      write(11,'(d20.12,1x,200(d20.12,1x))')
     &           t, (y(i),i=1,neq)
      ierr = tcgetms2cc( y, neq, ymfr ) ;
      write(12,'(d20.12,1x,200(d20.12,1x))')
     &           t, (ymfr(i),i=1,ns)
      ierr = tcgetms2hmixms( y, neq, ymfr(1) ) ;
      write(13,'(d20.12,1x,d20.12)') t, ymfr(1)

      close(10)
      close(11)
      close(12)
      close(13)
#endif
 
c-Output ignition delay
      open(unit=14,file='tid.dat')
      write(14,'(e20.12,1x,e20.12)')time_id1, time_id2
      close(14)

c-Reset the library
      call tcreset()

      if ( istate.lt.0)  goto 80
c      write(*,*)' Statistics : '
c      write(*,60)  iwork(11), iwork(12), iwork(13)
c   60 format(/' No. steps =',i6,',  No. f-s =',i6,',  No. J-s =',i6)

      stop
   80 write(*,90) istate
   90 format(///' Error halt.. istate =',I3)
      stop
      end

      include 'rhsjac.f'
c
c-Input file
c
      subroutine initrun(itermax, outfreq, mktab,deltatsta, deltatmax, 
     &                   tmax, Tempini, Temp_id, deltaTemp, smsfrin, 
     &                   snamein, nspecin, mechfile, thermofile)

      implicit none

      integer       itermax, outfreq, mktab, nspecin, ispec
      real*8        deltatsta, deltatmax, tmax, Tempini, Temp_id, 
     &              deltaTemp
      real*8        smsfrin(*)
      character*100 mechfile, thermofile
      character*18  snamein(*)

      character*11 keyname
      integer      itmp1
      real*8       rtmp1

      ispec = 0
 100  continue

        read(5,'(A11)',advance='no') keyname
        select case (keyname)
          case ("NiterMax")
            read(5,*)itmp1
            itermax = itmp1
            write(*,'(A11,A1,I7)')keyname,' ',itmp1
          case ("oFreq")
            read(5,*)itmp1
            outfreq = itmp1
            write(*,'(A11,A1,I7)')keyname,' ',itmp1
          case ("deltat")
            read(5,*)rtmp1
            deltatsta = rtmp1
            write(*,'(A11,A1,E14.6)')keyname,' ',rtmp1
          case ("deltatMax")
            read(5,*)rtmp1
            deltatmax = rtmp1
            write(*,'(A11,A1,E14.6)')keyname,' ',rtmp1
          case ("tEnd")
            read(5,*)rtmp1
            tmax = rtmp1
            write(*,'(A11,A1,E14.6)')keyname,' ',rtmp1
          case ("Tini")
            read(5,*)rtmp1
            Tempini = rtmp1
            write(*,'(A11,A1,E14.6)')keyname,' ',rtmp1
          case ("Temp_id")
            read(5,*)rtmp1
            Temp_id = rtmp1
            write(*,'(A11,A1,E14.6)')keyname,' ',rtmp1
          case ("deltaTemp")
            read(5,*)rtmp1
            deltaTemp = rtmp1
            write(*,'(A11,A1,E14.6)')keyname,' ',rtmp1
          case ("mech")
            read(5,*)mechfile
            write(*,'(A11,A1,A)')keyname,' ',mechfile
          case ("thermo")
            read(5,*)thermofile
            write(*,'(A11,A1,A)')keyname,' ',thermofile
          case ("withTab")
            read(5,*)itmp1
            mktab = itmp1
            write(*,'(A11,A1,I7)')keyname,' ',itmp1
          case ("spec")
            ispec = ispec+1
            read(5,*)snamein(ispec),smsfrin(ispec)
            write(*,'(A11,A1,A18,A1,E20.12)')keyname,' ',snamein(ispec),
     &            ' ',smsfrin(ispec)
          case ("END")
            write(*,'(A11)')keyname
          case default
            write(*,*)'No info on : ',keyname,' -> Abort'
            stop
        end select

      write(*,*) keyname
      if ( keyname.ne.'END') goto 100

      nspecin = ispec 

      return
      end
