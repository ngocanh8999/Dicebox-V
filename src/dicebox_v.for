c=======================================================================
c DICEBOX-V: A Versatile and High-Performance Multi-Threaded Extension
c            of the DICEBOX Statistical Model Code
c
c Description:
c   DICEBOX-V (where "V" stands for "Versatile") is an independent,
c   multi-threaded OpenMP fork of DICEBOX version 1.0.1 (originally
c   developed by Milan Krticka, Charles University, Prague, based on the
c   statistical model and algorithm by Frantisek Becvar, Nucl. Instrum.
c   Methods Phys. Res. A 417 (1998) 434).
c
c   This modified code is developed and maintained separately by the
c   DICEBOX-V authors (N. N. Anh et al.). It is unaffiliated with the
c   original DICEBOX authors, who neither maintain nor assume any
c   responsibility for this modified code.
c
c License:
c   This program is free software: you can redistribute it and/or modify
c   it under the terms of the GNU General Public License as published by
c   the Free Software Foundation, either version 3 of the License, or
c   (at your option) any later version.
c
c   This program is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c   GNU General Public License for more details.
c
c   You should have received a copy of the GNU General Public License
c   along with this program. If not, see <https://www.gnu.org/licenses/>.
c=======================================================================
c
      PROGRAM DICE_EVENT
c
      INTEGER*4    MAXIS,MBATCH
      PARAMETER    (MAXJC  = 49,   !maximum spin in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000,
     *              MBATCH = 2048            )
c
      LOGICAL        LOGS,ponverze,lopopgs,sidlev,ponverk,logfil,
     &               LPAIRED_GERMS,LNO_WCACHE
      INTEGER        NARGS, IARG
      CHARACTER*64   CARG
      COMMON /PGFLAG/ LPAIRED_GERMS
      COMMON /CWOPT/  LNO_WCACHE
c      CHARACTER*80   OUTNAME
      CHARACTER*12   NAME
      CHARACTER*3    EXT
      INTEGER*4      IR1,IR2,IR3,IR4,STEPS,pointers,IRCON,MULT,NTOTAL,
     &               IRCONc,ISCON,ISDIS,IRINIT,
     &               IR1_BASE,IR3_BASE,IR4_BASE,REAL_ID,SEED_EVENT,
     &               NTRUNC126
      INTEGER        DEKOD,DENUM,DELEV,DEPARITY
      REAL*8         TOTCON,GACON,STCON
      DIMENSION      TOTDIS(0:2),TOTCON(0:2),
     &               STCON(0:2,0:MAXBIN,-2:2,0:1),GACON(0:2,-2:2,0:1),
     &               GADIS(0:2,-2:2,0:1),STDIS(0:2,0:30,-2:2,0:1),
     &               ISCON(0:2,0:MAXBIN,-2:2,0:1),
     &               ISDIS(0:2,30,0:MAXJC,0:1)
      DIMENSION      POPULT(99),POPERT(99),COVAP(99,99),CONO(99,99),
     &               POPULS(99),POPERS(99),COVAS(99,99),OLDPOP(99)
      INTEGER*4      NSTEP_B(MBATCH)
      REAL*4         ELQQ_B(0:126,MBATCH),SPQQ_B(0:126,MBATCH),
     *               WIQQ_B(0:126,MBATCH),DMQQ_B(126,MBATCH)
      INTEGER*4      IPQQ_B(0:126,MBATCH),ICQQ_B(126,MBATCH)
      COMMON /EVBATCH/ ELQQ_B,SPQQ_B,IPQQ_B,ICQQ_B,DMQQ_B,WIQQ_B,
     *                 NSTEP_B
      PARAMETER (MAX_OMP_TH = 128)
      INTEGER*8      POPTLEV_INT(99),POPSLEV_INT(99),SPLITMIX64
      INTEGER*8      POPTLEV_TH(99,0:MAX_OMP_TH),
     &               POPSLEV_TH(99,0:MAX_OMP_TH)
      INTEGER        NDEAD_TH(0:MAX_OMP_TH),NTRUNC_TH(0:MAX_OMP_TH)
      INTEGER        TH_SYNC(0:MAX_OMP_TH)
      INTEGER        MY_ID,NTHREADS_ACT,ISYNC,ITH,NTH_MAX
      INTEGER        OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,
     &               OMP_GET_MAX_THREADS
      REAL*8         TOTCON_LOC,GACON_LOC,STCON_LOC
      DIMENSION      TOTDIS_LOC(0:2),TOTCON_LOC(0:2),
     &               STCON_LOC(0:2,0:MAXBIN,-2:2,0:1),
     &               GACON_LOC(0:2,-2:2,0:1),GADIS_LOC(0:2,-2:2,0:1),
     &               STDIS_LOC(0:2,0:30,-2:2,0:1),
     &               ISCON_LOC(0:2,0:MAXBIN,-2:2,0:1),
     &               ISDIS_LOC(0:2,30,0:MAXJC,0:1)
      CHARACTER*90   EV_LINE(100,MBATCH)
      INTEGER*4      EV_LEN(100,MBATCH),EV_NLIN(MBATCH)
      COMMON /EVSTRBATCH/ EV_LEN,EV_NLIN,EV_LINE
      INTEGER*4      NL_EV,J_START,J_END,N_CHUNK,I_L,J
      INTEGER*4      IBATCH,NBATCH,J_B,IR1_EV,IR3_EV,IR4_EV,
     &               NSTEP_EV,J_ST
c
C     Hereafter, article 'CON' is related to a 'CONtinuum' part of
C     the level system, while 'DIS' -- to the 'DIScrete' part.
C     Suffixes 'p' and 's' are related to primary and secondary
C     transitions, respectively.
C
      COMMON /GAU/   U,IFLAG
C$OMP THREADPRIVATE(/GAU/)
      COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd,dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
     &       /incr/  nfilev,STEPS,enrgfin(18),buff(251),
     &               ncum,cumwidth(14),edo(18),eup(18)
      COMMON /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
     &               NOPTCS,NLINc,CAPFR(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     &               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     &               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     &               PAR_E1(3),DIPSLP,DIPZER
     &       /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
      COMMON /kon/   ponverze,ponverk
C$OMP THREADPRIVATE(/kon/)
      COMMON /intr/  last,ir1,ir2,ir3,ir4,radwid,radwdi,
     &               fnoav,fnodi,rnn(18,14),
     &               rmult(18,126),rovfl(18),avea(18,14),erra(18,14),
     &               puvea(18,14),cova(18,18,14),aves(1:18,0:MAXBIN),
     &               errs(1:18,0:MAXBIN),
     &               diss(1:18,0:MAXBIN),poptlev(99),popslev(99)
     &       /readi/ nevents,NSREAL,NRL,ecrit,xrayk,xrayl,numlev,
     &               elowlev(99),elowsp(99),ilowip(99),
     &               STDISa(0:2,0:30,-2:2,0:1),RADW
     &       /grnd/  logs,kgs,LOpopGS,KpopGS
     &       /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
     &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
      COMMON /monit/ ELQQ(0:126),SPQQ(0:126),IPQQ(0:126),ICQQ(126),
     &               DMQQ(126),WIQQ(0:126)
C$OMP THREADPRIVATE(/monit/)
      COMMON /SFCE/  SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
     &               GE1DIS(2),GM1DIS(2),GE2DIS(2)
C$OMP THREADPRIVATE(/SFCE/)
      COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
     *               NEIGENVAL(0:MAXJC,0:1),LMODE
     &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS
     &       /rnd/   IRINIT(4,999)

c     lpreru=.FALSE.
      LOGS=.FALSE.
      lopopgs=.false.
      ponverze=.FALSE.
      ponverk=.FALSE.
      NDEAD=0
      NTRUNC126=0
      DO J=0,MAXJC
        DO K=0,1
          NDIS(J,K)=0
        ENDDO
      ENDDO
      DO K=0,1
        DO L=0,30
          DO M=-2,2
            DO MM=0,2
              STDIS(MM,L,M,K)=0.
              STDISa(MM,L,M,K)=0.
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO K=1,99
        DO L=0,20
          SALL(K,L)=0.0
        ENDDO
        POPULT(K)=0.0
        POPERT(K)=0.0
        POPULS(K)=0.0
        POPERS(K)=0.0
        DO L=1,99
          COVAP(K,L)=0.0
          COVAS(K,L)=0.0
        ENDDO
      ENDDO
      RADWID=0.0
      RADWDI=0.0
c

cc    ! initial reading
c      CALL get_command_argument(1, NAME)
c      IF (LEN_TRIM(NAME).EQ.0) THEN
c          STOP 'Invalid input: need an input name as the first argument'
c      ENDIF
c      NAME=TRIM(NAME)
c      INQUIRE(FILE=NAME, EXIST=logfil)
c      IF (.not.logfil) THEN
c        STOP 'Invalid input: file not found'
c      ENDIF
c      CALL READ_EV(NAME)

C
C     Parse command-line arguments (e.g. --paired-germs, --no-wcache)
C
      LPAIRED_GERMS = .FALSE.
      LNO_WCACHE = .FALSE.
      NARGS = COMMAND_ARGUMENT_COUNT()
      DO IARG = 1, NARGS
        CALL GET_COMMAND_ARGUMENT(IARG, CARG)
        IF (INDEX(CARG, '--paired-germs') .GT. 0) THEN
          LPAIRED_GERMS = .TRUE.
        ENDIF
        IF (INDEX(CARG, '--no-wcache') .GT. 0 .OR.
     *      INDEX(CARG, '--disable-wcache') .GT. 0) THEN
          LNO_WCACHE = .TRUE.
        ENDIF
      ENDDO
      IF (LNO_WCACHE) THEN
        WRITE(*,*)
     *    'INFO: Dynamic width cache (/WCACHE/) is DISABLED.'
        WRITE(*,*)
     *    '      All partial widths will be computed directly.'
      ENDIF

      NAME='DICE_EV.DAT'
      CALL READ_EV(NAME)
      CALL RANDINIT()
      IF (SPLITMIX64(0_8) .NE. -2152535657050944081_8) THEN
        WRITE(*,*) 'FATAL: SPLITMIX64 failed self-test!'
        WRITE(*,*)
     *    '64-bit integer overflow wrapping (-fwrapv) is broken.'
        STOP 'FATAL: SPLITMIX64 verification failed'
      ENDIF
      IF (SPLITMIX64(INT(NARGS, 8) + 123456789_8) .EQ. 0_8) THEN
        STOP 'FATAL: SPLITMIX64 dynamic check failed'
      ENDIF
      CALL WCACHE_INIT()
C
C     Check OpenMP thread count against MAX_OMP_TH
C
      NTH_MAX = 1
C$    NTH_MAX = OMP_GET_MAX_THREADS()
      IF (NTH_MAX .GT. MAX_OMP_TH) THEN
        WRITE(*,*) 'WARNING: Requested OpenMP threads (', NTH_MAX,
     *             ') exceeds maximum supported (', MAX_OMP_TH, ').'
        WRITE(*,*) '         Automatically clamping thread count ',
     *             'to ', MAX_OMP_TH
C$      CALL OMP_SET_NUM_THREADS(MAX_OMP_TH)
      ENDIF
C
C     Generation of GOE eigenvalues
c
      IF (LMODE.EQ.1) THEN
       NDIM=700              !Maximum allowed value is NDIM=1000
       NDIM=300              !Maximum allowed value is NDIM=1000
       CALL GENERATE_GOE_EIGEN_VAL(IR1,NDIM)
      ENDIF
c
c     PSFs used in tabular form
c
      OPEN (UNIT=12,FILE='PSF_GS.DAT',STATUS='UNKNOWN')
      OPEN (UNIT=13,FILE='PSF_INI.DAT',STATUS='UNKNOWN')

      EG_MAX = 11.0
      EG_STEP = 0.05
      NEG = INT (EG_MAX/EG_STEP)
      DO I = 1, NEG
        eg = FLOAT(I) * EG_STEP
        sf1_ini = sgamma(eg,bn,1)/eg**3     !from initial level
        sf2_ini = sgamma(eg,bn,3)/eg**3
        sf3_ini = sgamma(eg,bn,4)/eg**5
c        sf4_ini = sgamma(eg,bn,4)/eg**3
        sf1_fin = sgamma(eg,eg,1)/eg**3     !to GS
        sf2_fin = sgamma(eg,eg,3)/eg**3
        sf3_fin = sgamma(eg,eg,4)/eg**5
c        sf4_fin = sgamma(eg,eg,4)/eg**3
        write(12,201) eg,sf1_fin,sf2_fin,sf3_fin
        write(13,201) eg,sf1_ini,sf2_ini,sf3_ini
  201   format(f7.3,4e13.5)
      ENDDO
      CLOSE (12)
      CLOSE (13)
c
c          The loop through number of realisation of nucleus.
c
      DO NUC=1,NSREAL
c
c       write(*,*) NUC,IR1, IR2, IR3, IR4
        IR1 = IRINIT(1,NUC) !level scheme
        IR2 = IRINIT(2,NUC) !radiation widths in form of precursors and low-lying intensity fluctuations
        IR3 = IRINIT(3,NUC) !MC of cascades - actual search for the final state
        IR4 = IRINIT(4,NUC) !MC of cascades - the seed for a) mixing of primaries b) "coin flip" of mixing ratio
c
        CALL LEVELS (SPINC(1),IR1) ! Generation of level scheme
c
        IF (NTOTAL + nddd .GT. MAXIS) THEN
          WRITE (*,*)
     &      'The level density is too high. Increase IRCON dimension'
          WRITE (*,*) 'The required dimension is ', NTOTAL + nddd
          STOP 'ERROR'
        ENDIF
c
c       Writting levels 
c
        IF (ISWLS.EQ.1) THEN   ! Writing generated levels if requested
          WRITE (EXT,211) NUC
          S0=SPINC(1)-FLOAT(INT(SPINC(1)))
          IF (S0.LE.1e-3) THEN  ! only for security reason
           S0 = 0.0
          ELSE
           S0 = 0.5
          ENDIF
          OPEN (UNIT=13,FILE='LEVELS.'//EXT,STATUS='UNKNOWN')
          DO IS = 0, MAXJC
           SPFI = FLOAT(IS) + S0
           DO IPFI = 0, 1
            DO IBIN = 1, NBIN
              NL = LEVCON(IBIN,IS,IPFI)-LEVCON(IBIN-1,IS,IPFI)
              DO IL = 1, NL
                Q=1.-FLOAT(2*IL-1)/FLOAT(2*NL)
                ENRG = BN - (FLOAT(IBIN)-Q)*DELTA
                WRITE(13,212) ENRG, SPFI, IPFI
              ENDDO
            ENDDO ! IBIN
            NL = ndis(ISUBSC(spfi),ipfi)
            DO IL = 1, NL
              ENRG = endis(IL,ISUBSC(spfi),ipfi)
              WRITE(13,212) ENRG, SPFI, IPFI
            ENDDO
           ENDDO !IPFI
          ENDDO !IS
          CLOSE(13)
        ENDIF ! ISWLS
  211   FORMAT (I3.3)
  212   FORMAT(F8.5,F5.1,I2)  ! End of the part writing generated levels
c
        CALL GERMS(IR2) ! Precursors assigned to each level
        CALL WCACHE_CLEAR()
c
c       Intensities of low-lying transitions flucutuate
c
        NAME='DICE_EV.DAT'
        CALL READ_INT(NAME)

        DO IRL = 1, NRL
c
          CALL WCACHE_CLEAR()
          IF (ISWWR.EQ.1) CALL OPEN_IT(NUC,IRL)
          WRITE(*,*) ' Supra-Realisation : ',NUC,' (of',NSREAL,')'
          WRITE(*,*) ' Realisation : ',IRL,' (of',NRL,')'

C
C      The following DO loop serves for computing of mixing ratios
C      delta for primary transitions (E2 admixture is probably not
C      very important, so this is done in a very simple way - and
C      probably not fully correctly in the case NLINc=2)
C
          DO ilinc=1,nlinc
           DO ip=0,1
            DO is=0,MAXJC
             DO il=1,NDIS(is,ip)
              DO ipom=1,100
               dummy=ran0(IR4)
              ENDDO
              isdis(ilinc,il,is,ip)=IR4
             ENDDO
            ENDDO
           ENDDO
          ENDDO
c
c     Here, the procedure 'WIDTHS' is called in order to evaluate a total
c     radiative width and proper values of STDISp.
c
c     Label 19 - Security reason! If TOTDIS(one of capt. state)>1 
c     (possible due to the random distribution of intensities between c.s.)
c

   19	    CONTINUE
          IF (IPRIM.EQ.1) THEN
            NAME='DICE_EV.DAT'
            IF (NLINc.GT.1) CALL READ_AGAIN(NAME)
          ENDIF
 
          IREGI=0
          IBIN=0
          ILIN=1
          RADW=0.
c
c
          IF (IPRIM.EQ.1) THEN ! Known primaries
           DO ILINc=1,NLINc
            CALL WIDTHS (ILINc,IPINC,SPINC(ILINc),IBIN,ILIN,TOTCON,
     &                   STCON,GACON,ISCON,TOTDIS,STDISa,GADIS,ISDIS)
           ENDDO
           DO ILINc=1,NLINc
            IF(TOTDIS(ILINc).GE.1.) GOTO 19
           ENDDO	 

           DO ILINc=1, NLINc
            RADWp=sngl(TOTCON(ILINc))/(1-TOTDIS(ILINc))
            DO IPFI=0,1
             SP=SPINC(ILINc)-INT(SPINC(ILINc)+.25)
             DO IS = 0, MAXJC
              SPFI = FLOAT(IS) + SP
              ISBS=NINT(SPFI+.25)-NINT(SPINC(ILINc)+.25)
              IF ((ISBS.GE.-2).AND.(ISBS.LE.2)) THEN
               DO L=1,NDIS(ISUBSC(SPFI),IPFI)
                STDIS(ILINc,L,ISBS,IPFI)=STDISa(ILINc,L,ISBS,IPFI)*RADWp
               ENDDO
              ENDIF
             ENDDO ! IS / SPFI
            ENDDO ! IPFI
            RADW=RADW+RADWp*CAPFR(ILINc)
           ENDDO ! ILINc
           IREGI=0
           IBIN=0
           ILIN=1
C
           DO ILINc=1,NLINc
            TOTDIS(ILINc)=0.
            DO IS=-2,2
             DO IPFI=0,1
              GADIS(ILINc,IS,IPFI)=0.
             ENDDO
            ENDDO
            SP=SPINC(ILINc)-INT(SPINC(ILINc)+.25)
            ISP=INT(SPINC(ILINc)+.25)
            DO IPFI=0,1
             DO IS=ISP-2,ISP+2
              SPFI=SP+FLOAT(IS)
              ISBS=NINT(SPFI+.25)-NINT(SPINC(ILINc)+.25)
              IF (ITYPE(SPINC(ILINc),IPINC,SPFI,IPFI).NE.0) THEN
               TOTDIS(ILINc)=TOTDIS(ILINc)+
     &            STDIS(ILINc,NDIS(ISUBSC(SPFI),IPFI),ISBS,IPFI)
               GADIS(ILINc,ISBS,IPFI)=TOTDIS(ILINc)
              ENDIF
             ENDDO
            ENDDO
           ENDDO

          ELSE
           DO ILINc=1,NLINc   ! Uknown primary intensities (resonances)
            CALL WIDTHS(ILINc,IPINC,SPINC(ILINc),IBIN,ILIN,
     &                TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS)
           ENDDO

           DO ILINc=1, NLINc
             RADW=RADW+sngl(TOTCON(ILINc))+TOTDIS(ILINc)
           ENDDO
           IREGI=0
           IBIN=0
           ILIN=1
          ENDIF !IPRIM


c
c   Population of low-lying levels - setting initial (zero) values
c
         DO I=1,NUMLEV
           POPTLEV(I)=0.0
           POPSLEV(I)=0.0
         ENDDO
C
c**************************************
c       The master DO-loop
c
         WRITE(*,*) ' Event (N=',NEVENTS,') : '
         WRITE(*,*)
         IR1_BASE = IR1
         IR3_BASE = IR3
         IR4_BASE = IR4
         REAL_ID = IRL + (NUC - 1) * NRL
         DO I = 1, NUMLEV
           POPTLEV_INT(I) = 0_8
           POPSLEV_INT(I) = 0_8
         ENDDO

         DO IBATCH = 1, NEVENTS, MBATCH
           NBATCH = MIN(MBATCH, NEVENTS - IBATCH + 1)
           NTHREADS_ACT = 1
           DO ITH = 0, MAX_OMP_TH
             TH_SYNC(ITH) = 0
             DO I = 1, NUMLEV
               POPTLEV_TH(I, ITH) = 0_8
               POPSLEV_TH(I, ITH) = 0_8
             ENDDO
             NDEAD_TH(ITH) = 0
             NTRUNC_TH(ITH) = 0
           ENDDO

C$OMP PARALLEL
C$OMP& COPYIN(/MAINE/, /SFCE/)
C$OMP& SHARED(TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS,
C$OMP&        POPTLEV_TH,POPSLEV_TH,NDEAD_TH,NTRUNC_TH,TH_SYNC,
C$OMP&        NTHREADS_ACT)
C$OMP& PRIVATE(TOTCON_LOC,STCON_LOC,GACON_LOC,ISCON_LOC,
C$OMP&         TOTDIS_LOC,STDIS_LOC,GADIS_LOC,ISDIS_LOC,
C$OMP&         J_B,IEV,IR1_EV,IR3_EV,IR4_EV,sidlev,NSTEP_EV,pointers,
C$OMP&         ILINc,IBIN,ILIN,IPIN,SPIN,IPFI,SPFI,IBFI,ILFI,
C$OMP&         dmix2,sign,K,J_ST,NL_EV,J_START,J_END,N_CHUNK,J,MY_ID)
           MY_ID = 0
C$         MY_ID = OMP_GET_THREAD_NUM()
C$    IF (MY_ID .GT. MAX_OMP_TH) THEN
C$      WRITE(*,*) 'FATAL: Thread ID ', MY_ID, ' exceeds MAX_OMP_TH'
C$      STOP 'FATAL: OMP thread ID overflow'
C$    ENDIF
C$OMP MASTER
C$         NTHREADS_ACT = OMP_GET_NUM_THREADS()
C$OMP END MASTER

C$OMP DO SCHEDULE(DYNAMIC, 16)
           DO J_B = 1, NBATCH
             IEV = IBATCH + J_B - 1
             IR1_EV = SEED_EVENT(IR1_BASE, REAL_ID, IEV, 1)
             IR3_EV = SEED_EVENT(IR3_BASE, REAL_ID, IEV, 3)
             IR4_EV = SEED_EVENT(IR4_BASE, REAL_ID, IEV, 4)
             ponverze = .FALSE.
             ponverk = .FALSE.
             sidlev = .FALSE.
             EIN = BN
             NSTEP_EV = 1
             pointers = 0
             ILINc = 1
             IF (IVER.GE.1) THEN
               IF (NOPTCS.NE.1) THEN
                 IF (RAN0(IR1_EV).GT.CAPFR(1)) ILINc = 2
               ENDIF
             ENDIF
             ELQQ(0) = BN
             SPQQ(0) = SPINc(ILINc)
             IPQQ(0) = IPINc
             TOTCON_LOC(0) = 0.0D0
             TOTDIS_LOC(0) = 0.0
             IREGI = 0
             IBIN = 0
             ILIN = 1
             IF ((TOTCON(ILINc) + DBLE(TOTDIS(ILINc))).LE.0.0D0) THEN
               NDEAD_TH(MY_ID) = NDEAD_TH(MY_ID) + 1
               NSTEP_EV = 0
               DO J_ST = 1, 126
                 ELQQ(J_ST) = 0.0
                 SPQQ(J_ST) = 0.0
                 WIQQ(J_ST) = 0.0
                 DMQQ(J_ST) = 0.0
                 ICQQ(J_ST) = 0
               ENDDO
               GOTO 440
             ENDIF
             CALL ONESTEP(ILINc,IPINC,SPINC(ILINc),IBIN,ILIN,TOTCON,
     &                    STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS,
     &                    IPFI,SPFI,IBFI,ILFI,dmix2,sign,IR3_EV,IR4_EV)
             DO WHILE (EFI.GT.0.)
               pointers = pointers + 1
               ELQQ(NSTEP_EV) = EFI
               SPQQ(NSTEP_EV) = SPFI
               IPQQ(NSTEP_EV) = IPFI
               DMQQ(NSTEP_EV) = sign * sqrt(dmix2)
               IF (NSTEP_EV.EQ.1) THEN
                 WIQQ(NSTEP_EV) = TOTCON(ILINc) + TOTDIS(ILINc)
               ELSE
                 WIQQ(NSTEP_EV) = TOTCON_LOC(0) + TOTDIS_LOC(0)
               ENDIF
               IF (ELQQ(NSTEP_EV-1).LE.ECRIT) WIQQ(NSTEP_EV) = 0.0
               IF (IREGI.GT.0) THEN
                 K = DEKOD(ILFI, ISUBSC(SPFI), IPFI)
                 IF (K.GT.0 .AND. K.LE.NUMLEV) THEN
                   POPTLEV_TH(K, MY_ID) = POPTLEV_TH(K, MY_ID) + 1_8
                   IF (.NOT.SIDLEV) THEN
                     POPSLEV_TH(K, MY_ID) = POPSLEV_TH(K, MY_ID) + 1_8
                     SIDLEV = .TRUE.
                   ENDIF
                 ENDIF
               ENDIF
               IF (ponverze) THEN
                 IF (ponverk) THEN
                   ICQQ(NSTEP_EV) = 1
                 ELSE
                   ICQQ(NSTEP_EV) = 2
                 ENDIF
               ELSE
                 ICQQ(NSTEP_EV) = 0
               ENDIF
               ponverze = .FALSE.
               ponverk = .FALSE.
               IPIN = IPFI
               SPIN = SPFI
               IBIN = IBFI
               ILIN = ILFI
               IF (IREGI.LT.2) THEN
                 CALL WIDTHS(0,IPIN,SPIN,IBIN,ILIN,
     &              TOTCON_LOC,STCON_LOC,GACON_LOC,ISCON_LOC,
     &              TOTDIS_LOC,STDIS_LOC,GADIS_LOC,ISDIS_LOC)
                 IF ((SNGL(TOTCON_LOC(0))+TOTDIS_LOC(0)).LE.0.) THEN
                   NDEAD_TH(MY_ID) = NDEAD_TH(MY_ID) + 1
                   GOTO 440
                 ENDIF
               ENDIF
               IF (IREGI.EQ.2) THEN
                 IF (denum(dekod(ILFI,ISUBSC(SPFI),IPFI)).EQ.0 .OR.
     &               sall(dekod(ILFI,ISUBSC(SPFI),IPFI),
     &                    denum(dekod(ILFI,ISUBSC(SPFI),IPFI))).LE.0.0)
     &               THEN
                   NDEAD_TH(MY_ID) = NDEAD_TH(MY_ID) + 1
                   GOTO 440
                 ENDIF
               ENDIF
               IF (NSTEP_EV.GE.126) THEN
                 NTRUNC_TH(MY_ID) = NTRUNC_TH(MY_ID) + 1
                 GOTO 440
               ENDIF
               NSTEP_EV = NSTEP_EV + 1
               CALL ONESTEP(0,IPIN,SPIN,IBIN,ILIN,TOTCON_LOC,STCON_LOC,
     &           GACON_LOC,ISCON_LOC,TOTDIS_LOC,STDIS_LOC,GADIS_LOC,
     &           ISDIS_LOC,IPFI,SPFI,IBFI,ILFI,dmix2,sign,IR3_EV,IR4_EV)
             ENDDO
             ELQQ(NSTEP_EV) = EFI
             SPQQ(NSTEP_EV) = SPFI
             IPQQ(NSTEP_EV) = IPFI
             DMQQ(NSTEP_EV) = sign * sqrt(dmix2)
             IF (NSTEP_EV.EQ.1) THEN
               WIQQ(NSTEP_EV) = TOTCON(ILINc) + TOTDIS(ILINc)
             ELSE
               WIQQ(NSTEP_EV) = TOTCON_LOC(0) + TOTDIS_LOC(0)
             ENDIF
             IF (ELQQ(NSTEP_EV-1).LE.ECRIT) WIQQ(NSTEP_EV) = 0.0
             IF (lopopgs) THEN
               POPTLEV_TH(kpopgs, MY_ID) = POPTLEV_TH(kpopgs, MY_ID)+1_8
               IF (.NOT.sidlev) THEN
                 POPSLEV_TH(kpopgs, MY_ID) =
     &             POPSLEV_TH(kpopgs, MY_ID) + 1_8
               ENDIF
             ENDIF
             IF (ponverze) THEN
               IF (ponverk) THEN
                 ICQQ(NSTEP_EV) = 1
               ELSE
                 ICQQ(NSTEP_EV) = 2
               ENDIF
             ELSE
               ICQQ(NSTEP_EV) = 0
             ENDIF
             ponverze = .FALSE.
             ponverk = .FALSE.
  440        CONTINUE
             IF (ISWWR.EQ.1) THEN
               IF (ISWBN.EQ.0) THEN
                 NL_EV = 1
                 WRITE(EV_LINE(1,J_B),100) SPQQ(0),IPQQ(0),
     &                                     NSTEP_EV,ELQQ(0)
                 EV_LEN(1,J_B) = 24
                 IF (ISWEL.EQ.1) THEN
                   DO J_START = 1, NSTEP_EV, 8
                     J_END = MIN(NSTEP_EV, J_START + 7)
                     N_CHUNK = J_END - J_START + 1
                     NL_EV = NL_EV + 1
                     WRITE(EV_LINE(NL_EV,J_B),200)
     &                 (ELQQ(J), J = J_START, J_END)
                     EV_LEN(NL_EV,J_B) = 2 + N_CHUNK * 11
                   ENDDO
                 ENDIF
                 IF (ISWSP.EQ.1) THEN
                   DO J_START = 1, NSTEP_EV, 8
                     J_END = MIN(NSTEP_EV, J_START + 7)
                     N_CHUNK = J_END - J_START + 1
                     NL_EV = NL_EV + 1
                     WRITE(EV_LINE(NL_EV,J_B),300)
     &                 (SPQQ(J), J = J_START, J_END)
                     EV_LEN(NL_EV,J_B) = 2 + N_CHUNK * 11
                   ENDDO
                 ENDIF
                 IF (ISWPA.EQ.1) THEN
                   DO J_START = 1, NSTEP_EV, 8
                     J_END = MIN(NSTEP_EV, J_START + 7)
                     N_CHUNK = J_END - J_START + 1
                     NL_EV = NL_EV + 1
                     WRITE(EV_LINE(NL_EV,J_B),400)
     &                 (IPQQ(J), J = J_START, J_END)
                     EV_LEN(NL_EV,J_B) = 2 + N_CHUNK * 11
                   ENDDO
                 ENDIF
                 IF (ISWIC.EQ.1) THEN
                   DO J_START = 1, NSTEP_EV, 8
                     J_END = MIN(NSTEP_EV, J_START + 7)
                     N_CHUNK = J_END - J_START + 1
                     NL_EV = NL_EV + 1
                     WRITE(EV_LINE(NL_EV,J_B),400)
     &                 (ICQQ(J), J = J_START, J_END)
                     EV_LEN(NL_EV,J_B) = 2 + N_CHUNK * 11
                   ENDDO
                 ENDIF
                 IF (ISWMX.EQ.1) THEN
                   DO J_START = 1, NSTEP_EV, 8
                     J_END = MIN(NSTEP_EV, J_START + 7)
                     N_CHUNK = J_END - J_START + 1
                     NL_EV = NL_EV + 1
                     WRITE(EV_LINE(NL_EV,J_B),500)
     &                 (DMQQ(J), J = J_START, J_END)
                     EV_LEN(NL_EV,J_B) = 2 + N_CHUNK * 11
                   ENDDO
                 ENDIF
                 IF (ISWWI.EQ.1) THEN
                   DO J_START = 1, NSTEP_EV, 8
                     J_END = MIN(NSTEP_EV, J_START + 7)
                     N_CHUNK = J_END - J_START + 1
                     NL_EV = NL_EV + 1
                     WRITE(EV_LINE(NL_EV,J_B),500)
     &                 (WIQQ(J), J = J_START, J_END)
                     EV_LEN(NL_EV,J_B) = 2 + N_CHUNK * 11
                   ENDDO
                 ENDIF
                 EV_NLIN(J_B) = NL_EV
               ELSE
                 NSTEP_B(J_B) = NSTEP_EV
                 DO J_ST = 0, NSTEP_EV
                   ELQQ_B(J_ST, J_B) = ELQQ(J_ST)
                   SPQQ_B(J_ST, J_B) = SPQQ(J_ST)
                   IPQQ_B(J_ST, J_B) = IPQQ(J_ST)
                   WIQQ_B(J_ST, J_B) = WIQQ(J_ST)
                 ENDDO
                 DO J_ST = 1, NSTEP_EV
                   ICQQ_B(J_ST, J_B) = ICQQ(J_ST)
                   DMQQ_B(J_ST, J_B) = DMQQ(J_ST)
                 ENDDO
               ENDIF
             ENDIF
           ENDDO
C$OMP END DO
C$OMP ATOMIC WRITE RELEASE
           TH_SYNC(MY_ID) = 1
C$OMP END PARALLEL

           DO ITH = 0, NTHREADS_ACT - 1
C$OMP ATOMIC READ ACQUIRE
             ISYNC = TH_SYNC(ITH)
             DO I = 1, NUMLEV
               POPTLEV_INT(I) = POPTLEV_INT(I) + POPTLEV_TH(I, ITH)
               POPSLEV_INT(I) = POPSLEV_INT(I) + POPSLEV_TH(I, ITH)
             ENDDO
             NDEAD = NDEAD + NDEAD_TH(ITH)
             NTRUNC126 = NTRUNC126 + NTRUNC_TH(ITH)
           ENDDO

           IF (ISWWR.EQ.1) THEN
             IF (ISWBN.EQ.0) THEN
               DO J_B = 1, NBATCH
                 DO I_L = 1, EV_NLIN(J_B)
                   WRITE(12, '(A)') EV_LINE(I_L,J_B)(1:EV_LEN(I_L,J_B))
                 ENDDO
               ENDDO
             ELSE
               DO J_B = 1, NBATCH
                 STEPS = NSTEP_B(J_B)
                 DO J_ST = 0, STEPS
                   ELQQ(J_ST) = ELQQ_B(J_ST, J_B)
                   SPQQ(J_ST) = SPQQ_B(J_ST, J_B)
                   IPQQ(J_ST) = IPQQ_B(J_ST, J_B)
                   WIQQ(J_ST) = WIQQ_B(J_ST, J_B)
                 ENDDO
                 DO J_ST = 1, STEPS
                   ICQQ(J_ST) = ICQQ_B(J_ST, J_B)
                   DMQQ(J_ST) = DMQQ_B(J_ST, J_B)
                 ENDDO
                 CALL DO_IT(STEPS)
               ENDDO
             ENDIF
           ENDIF
         ENDDO

         IR1 = SEED_EVENT(IR1_BASE, REAL_ID, NEVENTS + 1, 1)
         IR3 = SEED_EVENT(IR3_BASE, REAL_ID, NEVENTS + 1, 3)
         IR4 = SEED_EVENT(IR4_BASE, REAL_ID, NEVENTS + 1, 4)
c*************************************
c
c*******    Cummulative quantities collected through the number
c                     of nuclear realizations :
c
         RNUC=FLOAT(IRL+(NUC-1)*NRL)
         RNUC1=RNUC-1.
         RNUC2=RNUC*RNUC
         RNUC12=RNUC1*RNUC1
c
c***    The total radiation width of the capturing state
c
         OLDWID=RADWID
         RADWID=(RNUC1*RADWID+RADW)/RNUC
         RADWDI=(RNUC1*(RADWDI+OLDWID*OLDWID)+RADW*RADW)/RNUC-
     &          RADWID*RADWID
c
c***    Population of low-lying levels (originaly nlowlev(level,nuc))
c
         DO K = 1, NUMLEV
           POPTLEV(K) = FLOAT(POPTLEV_INT(K)) / FLOAT(NEVENTS)
           POPSLEV(K) = FLOAT(POPSLEV_INT(K)) / FLOAT(NEVENTS)
           OLDPOP(K) = POPULT(K)
           POPULT(K) = (RNUC1*POPULT(K) + POPTLEV(K))/RNUC
           POPERT(K) = (RNUC1*(POPERT(K) + OLDPOP(K)*OLDPOP(K)) +
     *       POPTLEV(K)*POPTLEV(K))/RNUC - POPULT(K)*POPULT(K)
           DO L = 1, K
             COVAP(K,L) = (RNUC1*(COVAP(K,L) + OLDPOP(K)*OLDPOP(L)) +
     &         POPTLEV(K)*POPTLEV(L))/RNUC - POPULT(K)*POPULT(L)
           ENDDO
         ENDDO
         DO K = 1, NUMLEV
           OLDPOP(K) = POPULS(K)
           POPULS(K) = (RNUC1*POPULS(K) + POPSLEV(K))/RNUC
           POPERS(K) = (RNUC1*(POPERS(K) + OLDPOP(K)*OLDPOP(K)) +
     *       POPSLEV(K)*POPSLEV(K))/RNUC - POPULS(K)*POPULS(K)
           DO L = 1, K
             COVAS(K,L) = (RNUC1*(COVAS(K,L) + OLDPOP(K)*OLDPOP(L)) +
     &         POPSLEV(K)*POPSLEV(L))/RNUC - POPULS(K)*POPULS(L)
           ENDDO
         ENDDO
c
c***    Up-dating of the output files
c
         NAME='DICE.PRO'
         CALL WRITE_PARAMS(NAME)
         OPEN(9,FILE=NAME,ACCESS='APPEND',STATUS='OLD')
          WRITE(9,102) NOPTDE,NOPTE1,NOPTM1,NOPTE2
  102     FORMAT(' NOPTDE:',I2,'   NOPTE1:',I2,'   NOPTM1:',I2,
     &         '   NOPTE2:',I2)
          WRITE(9,*)
          WRITE(9,*) ' Real #',NUC,' of',NSREAL,
     &    '      Subr #',IRL,' of',NRL
          WRITE(9,*) ' Number of events',NEVENTS
          WRITE(9,*)
          ERR=SQRT(RADWDI)
          WRITE(9,*) 'Capt.state tot.rad.width (MeV): ',
     &              RADWID,' +/-',ERR
          WRITE (9,*)
          WRITE (9,*)
     &    'Number of cascades terminating at a dead end: ',NDEAD
          IF (NTRUNC126.GT.0) THEN
            WRITE(9,*)
     &      'WARNING: Cascades truncated at 126 steps:     ',NTRUNC126
            WRITE(*,*)
     &      'WARNING: Cascades truncated at 126 steps:     ',NTRUNC126
          ENDIF
c     ***  feeding of low-lying states ***
          WRITE(9,*) 'Population of low-lying states'
          if (rnuc1.GT.0) then   
	      rn1=rnuc/rnuc1
          else
	      rn1=1.
	    endif
          do i=1,numlev
           if (popert(i).GT.0.0) then
            poper=sqrt(popert(i)*rn1)
           else
            poper=0.0  
           endif
           WRITE(9,198)i,elowlev(i),popult(i),poper,
     &                elowsp(i),ilowip(i)
          enddo
  198     format(I3,f10.5,F11.5,' +- ',F7.5,F11.1,I3)
          write(9,*)
          WRITE(9,*)
     *        'Direct population of low-lying states from continuum'
          do i=1,numlev
           if (popers(i).GT.0.0) then
            poper=sqrt(popers(i)*rn1)
           else
            poper=0.0  
           endif
           WRITE(9,198)i,elowlev(i),populs(i),poper,
     &                 elowsp(i),ilowip(i)
          enddo
          write(9,*)
c
          WRITE(9,*) ' Population - covariance matrix'
          DO K=1,NUMLEV
           DO L=1,K
            CON=COVAP(K,K)*COVAP(L,L)
            IF (CON.GT.0.) THEN
              CONO(K,L)=COVAP(K,L)/SQRT(CON)
            ELSE
              CONO(K,L)=1.0
            ENDIF
           ENDDO
          ENDDO
          DO K=1,NUMLEV
c          WRITE(9,105) (CONO(K,L),L=1,NUMLEV)
       !    WRITE(9,*) (CONO(K,L),L=1,NUMLEV)   ! No covariance matrix written in the present version 
c  105     FORMAT(<NUMLEV>F7.3)
          ENDDO
          WRITE(9,*)

          WRITE(9,*) ' Sidefeeding - covariance matrix'
          DO K=1,NUMLEV
           DO L=1,K
            CON=COVAS(K,K)*COVAS(L,L)
            IF (CON.GT.0.) THEN
              CONO(K,L)=COVAS(K,L)/SQRT(CON)
            ELSE
              CONO(K,L)=1.0
            ENDIF
           ENDDO
          ENDDO
          DO K=1,NUMLEV
c          WRITE(9,105) (CONO(K,L),L=1,NUMLEV)
      !     WRITE(9,*) (CONO(K,L),L=1,NUMLEV)   ! No covariance matrix written in the present version 
          ENDDO

c        WRITE(9,*)
c        WRITE(9,*) 'last generator seed IR1: ',IR1
c        WRITE(9,*) '                    IR2: ',IR2
c        WRITE(9,*) '                    IR3: ',IR3
c        WRITE(9,*) '                    IR4: ',IR4

         CLOSE(9)
         CALL CLOSE_IT
c
c       write(*,*) NUC, IR1, IR2, IR3, IR4, 'l'

         IRCONc(1) = IR3
         IF (IR3 .LT. 2146483644) THEN
           IRCONc(2) = IR3 + 1000003
         ELSE
           IRCONc(2) = IR3 - 1000003
         ENDIF
         IF (IRCONc(2) .LE. 0) IRCONc(2) = 123456789

        ENDDO ! IRL
c   55 CONTINUE
      ENDDO ! NUC
c
c******************************************************************
  100 FORMAT (F5.1,I3,I5,F11.5)
  200 FORMAT (2X,8F11.5)
  300 FORMAT (2X,8F11.1)
  400 FORMAT (2X,8I11)
  500 FORMAT (2X,8E11.3)
  600 FORMAT (F5.1,I3,I5,127F11.5)
      CALL WCACHE_DESTROY()
c     STOP 'Execution of DICEBOX ended O.K.'
      STOP 'Execution of DICEBOX-V ended O.K.'
      END
c
C**********************************************************************
      FUNCTION ALPH_TOT (EI,SPI,IPI,EF,SPF,IPF,DMISQ,ITAB)
C**********************************************************************
C
C     Evaluates total or K-shell internal conversion coefficient using
C     precomputed log-energy grids and 4-point Lagrange weights (AICCX)
C     for ITAB = 1 (total ICC) or ITAB = 2 (K-shell ICC).
C
C           DMISQ   - SQUARED (!) mixing amplitude
C
C           MAEL0   - the type of the lowest order of gamma radiation.
C                     for MAgnetic radiation MAEL0=1, for ELectric
C                     MAEL0=0
C
C           MUL0    - the lowest multipolarity contributing to the
C                     transition (say, 2 in case of E2+M3)
C
C           MUL1    - the next contributing multipolarity
C                     (3 in this case)
C
      MUL0=NINT(ABS(SPI-SPF))
      IF (MUL0.EQ.0) MUL0=1
      MUL1=MUL0+1
      IAUX=(-1)**(IPI+IPF+MUL0+1)
      IF (IAUX.EQ.-1) THEN
         MAEL0=0    ! The dominating type is electric
         ELSE
         MAEL0=1    ! ... magnetic
      ENDIF
C
C     MAEL1 - the type  of the next-order contributing radiation
C
      MAEL1=1-MAEL0
      EG=EI-EF
      CT0=AICCX(EG,ITAB,MAEL0,MUL0)
      CT1=AICCX(EG,ITAB,MAEL1,MUL1)
c
      ALPH_TOT=(CT0+CT1*DMISQ)/(1.+DMISQ)
      RETURN
      END
C
      SUBROUTINE ALPH_PREP(EI,SPI,IPI,EF,SPF,IPF,ITAB,CT0,CT1)
C**********************************************************************
C     Precomputes loop-invariant ICC factors CT0 and CT1 for a bin.
      MUL0=NINT(ABS(SPI-SPF))
      IF (MUL0.EQ.0) MUL0=1
      MUL1=MUL0+1
      IAUX=(-1)**(IPI+IPF+MUL0+1)
      IF (IAUX.EQ.-1) THEN
         MAEL0=0
      ELSE
         MAEL0=1
      ENDIF
      MAEL1=1-MAEL0
      EG=EI-EF
      CT0=AICCX(EG,ITAB,MAEL0,MUL0)
      CT1=AICCX(EG,ITAB,MAEL1,MUL1)
      RETURN
      END
C
C***********************************************************************
      FUNCTION AICC(ETRA,TABEN,TABICC,MAEL,MUL,N)
C***********************************************************************
C
C     This subroutine provides cubic interpolation of values
C     Y=TABICC(MAEL,MUL,I) for fixed MAEL and MUL and variable I.
C     These values Y for various I are assumed to represent internal
C     conversion coefficients (ICC's) for generally non-equdistant
C     transition energies that are specified by TABEN(I). MAEL stands for
C     the type of a transition (0 for M1 and 1 for E1) while MUL means
C     multipolarity. N is the length of the table for a fixed MAEL and MUL.
C     Be careful, if ETRA (the energy at which it is desirable to get
C     the ICC-coefficient) is lower than TABEN(1) or higher than TABEN(N),
C     then the interpolation changes to extrapolation and difficulties may
C     start ...
C                                                    Version from 6-OCT-95
C
C***********************************************************************
C
      DIMENSION TABEN(50),TABICC(0:1,5,50),XX(4),YY(4),A(4)
C
      AICC=0.0
      IF (ETRA.LT.TABEN(1)) THEN
        AICC = TABICC(MAEL,MUL,1)
        GOTO 2
      ENDIF             ! ETRA=TABEN(1)
      IF (ETRA.GT.TABEN(N)) ETRA=TABEN(N)
      IF (ETRA.LE.TABEN(2)) THEN
         K=0
         ELSE
         IF (ETRA.LE.TABEN(N-1)) THEN
            ILO=3
            IHI=N-1
            DO WHILE (ILO.LT.IHI)
               IMID=(ILO+IHI)/2
               IF (ETRA.LE.TABEN(IMID)) THEN
                  IHI=IMID
               ELSE
                  ILO=IMID+1
               ENDIF
            ENDDO
            K=ILO-3
            GO TO 1
            ELSE
            K=N-4
         ENDIF
      ENDIF
C
    1 etrapom=etra
      etra=log(etra)
      DO J=1,4
         XX(J)=log(TABEN(K+J))
         YY(J)=log(TABICC(MAEL,MUL,K+J))
      ENDDO !J
      AICC=0.
      DO I=1,4
         A(I)=YY(I)
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)/(XX(J)-XX(I))
         ENDDO !J
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)*(XX(J)-ETRA)
         ENDDO !J
         AICC=AICC+A(I)
      ENDDO !I
C
      aicc=exp(aicc)
      etra=etrapom
    2 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE AICC_PRECOMP()
C***********************************************************************
C     Precomputes log(energy) and 4-point cubic Lagrange weights
C     A0(I) = log(Y_I) / prod_{J/=I}(X_J - X_I) for both ICC tables
C     (ITAB=1: total, ITAB=2: K-shell) in the exact division order of
C     AICC, preserving bit-identical interpolation values.
C
      COMMON /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
     &               CONVk(0:1,5,50)
      COMMON /AICCPC/ AIC_EN(50,2),AIC_XL(50,2),AIC_Y1(0:1,5,2),
     *                AIC_CF(4,0:46,0:1,5,2),NAIC(2)
      DIMENSION XX(4),YY(4)
      DO ITAB = 1, 2
        IF (ITAB.EQ.1) THEN
          N = NENT
        ELSE
          N = NENK
        ENDIF
        NAIC(ITAB) = N
        DO K = 1, MIN(N,50)
          IF (ITAB.EQ.1) THEN
            AIC_EN(K,ITAB) = ELENt(K)
          ELSE
            AIC_EN(K,ITAB) = ELENk(K)
          ENDIF
          AIC_XL(K,ITAB) = log(AIC_EN(K,ITAB))
        ENDDO
        DO MAEL = 0, 1
          DO MUL = 1, 5
            IF (ITAB.EQ.1) THEN
              AIC_Y1(MAEL,MUL,ITAB) = CONVt(MAEL,MUL,1)
            ELSE
              AIC_Y1(MAEL,MUL,ITAB) = CONVk(MAEL,MUL,1)
            ENDIF
          ENDDO
        ENDDO
        IF (N.LT.4 .OR. N.GT.50) GOTO 10
        DO K = 0, N-4
          DO MAEL = 0, 1
            DO MUL = 1, 5
              DO J = 1, 4
                XX(J) = AIC_XL(K+J,ITAB)
                IF (ITAB.EQ.1) THEN
                  YY(J) = log(CONVt(MAEL,MUL,K+J))
                ELSE
                  YY(J) = log(CONVk(MAEL,MUL,K+J))
                ENDIF
              ENDDO
              DO I = 1, 4
                A = YY(I)
                DO J = 1, 4
                  IF (I.NE.J) A = A/(XX(J)-XX(I))
                ENDDO
                AIC_CF(I,K,MAEL,MUL,ITAB) = A
              ENDDO
            ENDDO
          ENDDO
        ENDDO
   10   CONTINUE
      ENDDO
      RETURN
      END
C
C***********************************************************************
      FUNCTION AICCX(ETRA,ITAB,MAEL,MUL)
C***********************************************************************
C     Evaluates cubic log-Lagrange ICC interpolation from precomputed
C     weights in COMMON /AICCPC/, reproducing AICC bit-for-bit.
C
      COMMON /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
     &               CONVk(0:1,5,50)
      COMMON /AICCPC/ AIC_EN(50,2),AIC_XL(50,2),AIC_Y1(0:1,5,2),
     *                AIC_CF(4,0:46,0:1,5,2),NAIC(2)
      DIMENSION XX(4)
C
      N = NAIC(ITAB)
      IF (N.LT.4 .OR. N.GT.50) THEN
        IF (ITAB.EQ.1) THEN
          AICCX = AICC(ETRA,ELENt,CONVt,MAEL,MUL,NENT)
        ELSE
          AICCX = AICC(ETRA,ELENk,CONVk,MAEL,MUL,NENK)
        ENDIF
        RETURN
      ENDIF
      IF (ETRA.LT.AIC_EN(1,ITAB)) THEN
        AICCX = AIC_Y1(MAEL,MUL,ITAB)
        RETURN
      ENDIF
      IF (ETRA.GT.AIC_EN(N,ITAB)) ETRA=AIC_EN(N,ITAB)
      IF (ETRA.LE.AIC_EN(2,ITAB)) THEN
        K=0
      ELSEIF (ETRA.LE.AIC_EN(N-1,ITAB)) THEN
        ILO=3
        IHI=N-1
        DO WHILE (ILO.LT.IHI)
          IMID=(ILO+IHI)/2
          IF (ETRA.LE.AIC_EN(IMID,ITAB)) THEN
            IHI=IMID
          ELSE
            ILO=IMID+1
          ENDIF
        ENDDO
        K=ILO-3
      ELSE
        K=N-4
      ENDIF
      ELOG = log(ETRA)
      DO J=1,4
        XX(J) = AIC_XL(K+J,ITAB)
      ENDDO
      SUM = 0.
      DO I=1,4
        A = AIC_CF(I,K,MAEL,MUL,ITAB)
        DO J=1,4
          IF (I.NE.J) A = A*(XX(J)-ELOG)
        ENDDO
        SUM = SUM + A
      ENDDO
      AICCX = exp(SUM)
      RETURN
      END FUNCTION
C
c This file contain basic functions used in various programs concerning
c DICEBOX
c Namely: ISUBSC (SPIN)
c         ITYPE (SPIN,IPIN,SPFI,IPFI)
c         NPOISS (IR,AM)
c         GAUSS (IR)
c
C***********************************************************************
      FUNCTION ISUBSC(A)
C
C     The function ISUBSC is needed in such situations
C     when values of spin are to be transformed to values
C     that can serve as subscripts. The function converts
C     values of a real variable A to integer values
C     according to the following rule:
C
C                A = 0.0  ISUBSC = 0
C                    0.5           0
C                    1.0           1
C                    1.5           1
C                    2.0           2
C                    2.5           2
C                     .            .
C                     .            .
C                    etc.         etc.
C
C     This function is to be immune against a finite precission in
C     manipulations with real values.
C
C***********************************************************************
c
      IF (A.GE.0.0) THEN
        ISUBSC=INT(A+0.25)
      ELSE
        I=INT(A)
        X=A-FLOAT(I)
        IF (ABS(X-.5)-0.25.LT.0.0) THEN
C         The case of half-integer spin:
          ISUBSC=I
        ELSE
C         The case of integer spin:
          ISUBSC=NINT(A)
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
C***********************************************************************
      FUNCTION ITYPE(SPIN,IPIN,SPFI,IPFI)
C
C     The function ITYPE determines the type of gamma-transition.
C     The meaning of the variables is evident.
C      ITYPE=1: Pure E1-transition
C            2: Mixed (M1+E2)-transition
C            3: Pure M1-transition
C            4: Pure E2-transition
C     This function also checks wheather the spin of the final state
C     falls within the limits 0. to 50. or -- in case of odd product
C     nuclei -- within the limits 0.5 to 50.5.
C
C***********************************************************************
C
      JL=NINT(ABS(2.*SPIN-2.*SPFI))
      JU=NINT(2.*SPIN+2.*SPFI)
      IF((SPFI.LT.(-.25)).OR.(SPFI.GT.(50.75))) GOTO 3
      IF ((JL.LE.2).AND.(JU.GE.2)) GOTO 1
      IF ((JL.LE.4).AND.(JU.GE.4)) GOTO 2
    3 ITYPE=0
      RETURN
    2 IF (IPIN.NE.IPFI) GOTO 3
      ITYPE=4
      RETURN
    1 IF ((JL.LE.4).AND.(JU.GE.4)) GOTO 4
      IF (IPIN.EQ.IPFI) GOTO 5
    6 ITYPE=1
      RETURN
    5 ITYPE=3
      RETURN
    4 IF (IPIN.NE.IPFI) GOTO 6
      ITYPE=2
      RETURN
      END
C
C
C
C***********************************************************************
      FUNCTION NPOISS(IR,AM)
C
C     Poisson's random numbers
c
C***********************************************************************
C
      INTEGER*4 IR
      COMMON /GAU/ U,IFLAG
C$OMP THREADPRIVATE(/GAU/)
      IF (AM.GT.50.) GOTO 1
      Q=1.
      NPOISS=0
    2 Q=Q*RAN0(IR)
      IF (Q.LT.EXP(-AM)) RETURN
      NPOISS=NPOISS+1
      GOTO 2
    1 NPOISS=INT(AM+SQRT(AM)*GAUSS(IR))
      RETURN
      END
C
C
C
C***********************************************************************
      FUNCTION GAUSS(IR)
C
C     Normally distributed random numbers with
C     a zero mean and a unit variance
c
C***********************************************************************
C
      INTEGER*4 IR
      COMMON /GAU/ U,IFLAG
C$OMP THREADPRIVATE(/GAU/)
      IF(IFLAG) 2,1,2
    1 X1=2.*RAN0(IR)-1.
      X2=2.*RAN0(IR)-1.
      S=X1*X1+X2*X2
      IF (S.GE.1.) GOTO 1
      S=SQRT(-2.*ALOG(S)/S)
      U=X1*S
      IFLAG=1
      GAUSS=X2*S
      RETURN
    2 IFLAG=0
      GAUSS=U
      RETURN
      END
C
C***********************************************************************
      INTEGER*8 FUNCTION SPLITMIX64(X)
C
C     SplitMix64 64-bit pseudo-random generator / bit mixer.
C     Standard constants:
C       GAMMA = 0x9E3779B97F4A7C15 = -7046029254386353131_8
C       MUL1  = 0xBF58476D1CE4E5B9 = -4658895280553007687_8
C       MUL2  = 0x94D049BB133111EB = -7723592293110705685_8
C
C***********************************************************************
      INTEGER*8 X, Z
      INTEGER*8 C_GAMMA, C_MUL1, C_MUL2
      PARAMETER (C_GAMMA = -7046029254386353131_8)
      PARAMETER (C_MUL1  = -4658895280553007687_8)
      PARAMETER (C_MUL2  = -7723592293110705685_8)
C
      Z = X + C_GAMMA
      Z = IEOR(Z, ISHFT(Z, -30)) * C_MUL1
      Z = IEOR(Z, ISHFT(Z, -27)) * C_MUL2
      Z = IEOR(Z, ISHFT(Z, -31))
      SPLITMIX64 = Z
      RETURN
      END
C
C***********************************************************************
      INTEGER*4 FUNCTION SEED_EVENT(IBASE, REAL_ID, IEV, ISTRM)
C
C     Generates deterministic 31-bit seed in [1, 2147483646]
C     for event IEV of realization REAL_ID from base seed IBASE
C     and stream ID ISTRM.
C     Uses two rounds of SplitMix64 mixing to guarantee full avalanche
C     and zero cross-correlation between streams, events and realizations.
C
C***********************************************************************
      INTEGER*4 IBASE, REAL_ID, IEV, ISTRM
      INTEGER*8 K, Z, R64
      INTEGER*8 SPLITMIX64
      INTEGER*8 C_GAMMA, C_DOMAIN
      PARAMETER (C_GAMMA  = -7046029254386353131_8)
      PARAMETER (C_DOMAIN = 6364136223846793005_8)
C
C     Round 1: Combine base seed, realization ID, and stream ID
      K = IEOR(INT(IBASE, 8), ISHFT(INT(REAL_ID, 8), 32))
      K = K + INT(ISTRM, 8) * C_DOMAIN
      Z = SPLITMIX64(K)
C
C     Round 2: Mix event index IEV with golden gamma
      Z = Z + INT(IEV, 8) * C_GAMMA
      R64 = SPLITMIX64(Z)
C
C     Map to valid Park-Miller seed in [1, 2147483646]
      SEED_EVENT = INT(MOD(IAND(R64, 9223372036854775807_8),
     &                     2147483646_8)) + 1
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE GAUSS_SQ_PAIR(IR, GSQ1, GSQ2)
C
C     Generates a pair of independent chi-squared variates (1 d.o.f.)
C     using Marsaglia polar method without SQRT and without /GAU/.
C     Evaluates Park-Miller steps via 64-bit Mersenne folding and the
C     two-step jump multiplier 16807**2 mod (2**31-1) = 282475249.
C
C***********************************************************************
      INTEGER*4 IR
      INTEGER*8 X, P, U1, U2, M8, A1, A2
      PARAMETER (M8=2147483647_8, A1=16807_8, A2=282475249_8)
      REAL*4 GSQ1, GSQ2, X1, X2, S, FAC
      REAL*4 AM
      PARAMETER (AM = 2.0 / 2147483647.0)
C
      X = INT(IR, 8)
    1 P = A1 * X
      U1 = IAND(P, M8) + ISHFT(P, -31)
      IF (U1 .GE. M8) U1 = U1 - M8
      P = A2 * X
      U2 = IAND(P, M8) + ISHFT(P, -31)
      IF (U2 .GE. M8) U2 = U2 - M8
      X = U2
      X1 = REAL(INT(U1, 4)) * AM - 1.0
      X2 = REAL(INT(U2, 4)) * AM - 1.0
C
      S = X1*X1 + X2*X2
      IF (S .GE. 1.0 .OR. S .LE. 0.0) GOTO 1
      IR = INT(X, 4)
      FAC = -2.0 * ALOG(S) / S
      GSQ1 = X1 * X1 * FAC
      GSQ2 = X2 * X2 * FAC
      RETURN
      END
C
C***********************************************************************
      FUNCTION GAUSS_SQ_ONE(IR)
C
C     Generates a single chi-squared variate (1 d.o.f.)
C     using Marsaglia polar method without SQRT and without /GAU/,
C     with 64-bit Mersenne-folded jump-ahead Park-Miller updates.
C
C***********************************************************************
      INTEGER*4 IR
      INTEGER*8 X, P, U1, U2, M8, A1, A2
      PARAMETER (M8=2147483647_8, A1=16807_8, A2=282475249_8)
      REAL*4 GAUSS_SQ_ONE, X1, X2, S, FAC
      REAL*4 AM
      PARAMETER (AM = 2.0 / 2147483647.0)
C
      X = INT(IR, 8)
    1 P = A1 * X
      U1 = IAND(P, M8) + ISHFT(P, -31)
      IF (U1 .GE. M8) U1 = U1 - M8
      P = A2 * X
      U2 = IAND(P, M8) + ISHFT(P, -31)
      IF (U2 .GE. M8) U2 = U2 - M8
      X = U2
      X1 = REAL(INT(U1, 4)) * AM - 1.0
      X2 = REAL(INT(U2, 4)) * AM - 1.0
C
      S = X1*X1 + X2*X2
      IF (S .GE. 1.0 .OR. S .LE. 0.0) GOTO 1
      IR = INT(X, 4)
      FAC = -2.0 * ALOG(S) / S
      GAUSS_SQ_ONE = X1 * X1 * FAC
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE CHI2_SUM1(IR, NL, DSUM)
C
C     Fused register-resident accumulation of NL independent chi2_1
C     variates, bit-identical to summing DBLE(GAUSS_SQ_ONE(IR)).
C
C***********************************************************************
      INTEGER*4 IR, NL, IL
      REAL*8 DSUM, D
      INTEGER*8 X, P, U1, U2, M8, A1, A2
      PARAMETER (M8=2147483647_8, A1=16807_8, A2=282475249_8)
      REAL*4 X1, X2, S, FAC, AM
      PARAMETER (AM = 2.0 / 2147483647.0)
C
      IF (NL .LE. 0) RETURN
      X = INT(IR, 8)
      D = DSUM
      DO IL = 1, NL
    1   P = A1 * X
        U1 = IAND(P, M8) + ISHFT(P, -31)
        IF (U1 .GE. M8) U1 = U1 - M8
        P = A2 * X
        U2 = IAND(P, M8) + ISHFT(P, -31)
        IF (U2 .GE. M8) U2 = U2 - M8
        X = U2
        X1 = REAL(INT(U1, 4)) * AM - 1.0
        X2 = REAL(INT(U2, 4)) * AM - 1.0
        S = X1*X1 + X2*X2
        IF (S .GE. 1.0 .OR. S .LE. 0.0) GOTO 1
        FAC = -2.0 * ALOG(S) / S
        D = D + DBLE(X1 * X1 * FAC)
      ENDDO
      DSUM = D
      IR = INT(X, 4)
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE CHI2_SUM2(IR, NL, DSUM1, DSUM2)
C
C     Fused register-resident accumulation of NL paired chi2_1 variates,
C     bit-identical to calling GAUSS_SQ_PAIR(IR, G1, G2) in a loop.
C
C***********************************************************************
      INTEGER*4 IR, NL, IL
      REAL*8 DSUM1, DSUM2, D1, D2
      INTEGER*8 X, P, U1, U2, M8, A1, A2
      PARAMETER (M8=2147483647_8, A1=16807_8, A2=282475249_8)
      REAL*4 X1, X2, S, FAC, AM
      PARAMETER (AM = 2.0 / 2147483647.0)
C
      IF (NL .LE. 0) RETURN
      X = INT(IR, 8)
      D1 = DSUM1
      D2 = DSUM2
      DO IL = 1, NL
    1   P = A1 * X
        U1 = IAND(P, M8) + ISHFT(P, -31)
        IF (U1 .GE. M8) U1 = U1 - M8
        P = A2 * X
        U2 = IAND(P, M8) + ISHFT(P, -31)
        IF (U2 .GE. M8) U2 = U2 - M8
        X = U2
        X1 = REAL(INT(U1, 4)) * AM - 1.0
        X2 = REAL(INT(U2, 4)) * AM - 1.0
        S = X1*X1 + X2*X2
        IF (S .GE. 1.0 .OR. S .LE. 0.0) GOTO 1
        FAC = -2.0 * ALOG(S) / S
        D1 = D1 + DBLE(X1 * X1 * FAC)
        D2 = D2 + DBLE(X2 * X2 * FAC)
      ENDDO
      DSUM1 = D1
      DSUM2 = D2
      IR = INT(X, 4)
      RETURN
      END
C
C
C***********************************************************************
c     version 2.61                                         19.1.2005              
      FUNCTION DENSITY(EEXC,SPIN,IPAR)
C
C     Explicit expressions for level density
C     NOPTDE=0: CTF-model
C           =1: Bethe's level-density formula following formulation
C               of T.von Egidy et al., Nucl.Phys. A (1988)
c           =2: modified BSFG
c           =3: modified CTF
c         =4,5: BSFG with modified spin cut-off parameter
C
c     Changed a factor 1/2 in the row DENSITY=DENSITY*FJ*.5
c          => the parity-dependent level density allowed (PRC67,015803)
c
C***********************************************************************
C
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
     *       /denp/  LDENP,ZNUM,DENPC,DENPEN(3),IPDEN
C
      DENSITY=0.
      SIGSQ=0.0
      PAIRS=0.0
      IF (NOPTDE.EQ.0) THEN                                 ! CTF
        EEFF=EEXC-EZERO
        IF (EEFF.LE.0.) RETURN
        DENSITY=EXP(EEFF/TEMPER)/TEMPER
        SIGSQ=(.98*AMASS**.29)**2.
c        SIGSQ=(2.*AMASS**.29)**2.
      ELSEIF (NOPTDE.EQ.1) THEN                             ! BSFG
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0888*SQRT(ASHELL*EEFF)*AMASS**.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)
      ELSEIF (NOPTDE.EQ.2) THEN                             ! modified BSFG
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0888*SQRT(ASHELL*EEFF)*AMASS**.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)
        IF ((EEXC.GE.DENLO).AND.(EEXC.LE.DENHI)) THEN
         FCTDEN=DENPAR(1)+DENPAR(2)*EEXC+
     +          DENPAR(3)*EEXC**2+DENPAR(4)*EEXC**3
         DENSITY=DENSITY*FCTDEN
        ENDIF
      ELSEIF (NOPTDE.EQ.3) THEN                              ! modified CTF
        EEFF=EEXC-EZERO
        IF (EEFF.LE.0.) RETURN
        DENSITY=EXP(EEFF/TEMPER)/TEMPER
        SIGSQ=(.98*AMASS**.29)**2.
        IF ((EEXC.GE.DENLO).AND.(EEXC.LE.DENHI)) THEN
         FCTDEN=DENPAR(1)+DENPAR(2)*EEXC+
     +          DENPAR(3)*EEXC**2+DENPAR(4)*EEXC**3
         DENSITY=DENSITY*FCTDEN
        ENDIF
      ELSEIF (NOPTDE.EQ.4) THEN                             ! BSFG, s=0.1446 (Paar,Al-Quraishi)
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.1446*SQRT(ASHELL*EEFF)*AMASS**.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)
      ELSEIF (NOPTDE.EQ.5) THEN                             ! BSFG, another s (Al-Quraishi)
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0145*0.8*SQRT(EEFF/ASHELL)*AMASS**1.66667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)

      ELSEIF (NOPTDE.EQ.6) THEN                             ! BSFG - Von Egidy cut-off
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0146*(1+SQRT(1+4*ASHELL*EEFF))/2./ASHELL*AMASS**1.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)

      ELSEIF (NOPTDE.EQ.8) THEN                                 ! CTF - Von Egidy 2009
        EEFF=EEXC-EZERO09
        IF (EEFF.LE.0.) RETURN
        IF (EEXC.LT.0.5*PAIRING09) RETURN
        SIGSQ=.391*AMASS**0.675*(EEXC-0.5*PAIRING09)**0.312 
        DENSITY=EXP(EEFF/TEMPER09)/TEMPER09
      ELSEIF (NOPTDE.EQ.9) THEN                             ! BSFG - Von Egidy 2009
        EEFF=EEXC-EONE09
        IF (EEFF.LE.0.) RETURN
        IF (EEXC.LT.0.5*PAIRING09) RETURN
        SIGSQ=.391*AMASS**0.675*(EEXC-0.5*PAIRING09)**0.312 
        DENSITY=EXP(2.*SQRT(ASHELL09*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL09**.25*EEFF**1.25)

      ELSEIF (NOPTDE.EQ.11) THEN                            ! Goriely
        IF (EEXC.LE.0.) RETURN
        DENSITY = ALD(EEXC,SPIN,IPAR)
        GOTO 31
      ELSEIF (NOPTDE.EQ.12) THEN                            ! Kawano
        IF (EEXC.LE.0.) RETURN
        DENSITY = ALD(EEXC,SPIN,IPAR)
        GOTO 31
      ELSEIF (NOPTDE.EQ.14) THEN                            ! Goriely with T
        IF (EEXC.LE.0.) RETURN
        DENSITY = ALD(EEXC,SPIN,IPAR)
        GOTO 31

      ELSEIF (NOPTDE.EQ.7) THEN                             ! Voinov Mo BSFG
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0146*SQRT(EEFF/ASHELL)*AMASS**1.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)

      ELSEIF (NOPTDE.EQ.13) THEN                             ! BSFG - Von Egidy cut-off
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0146*(1+SQRT(1+4*ASHELL*EEFF))/2./ASHELL*AMASS**1.666667
        EEFF=EEXC-EZERO
        DENSITY=EXP(EEFF/TEMPER)/TEMPER

      ENDIF
c
      FJ_LOC=(SPIN+.5)*EXP(-(SPIN+.5)**2/(2.*SIGSQ))/SIGSQ
      IF ((NOPTDE.EQ.8).OR.(NOPTDE.EQ.9)) THEN
        IF (MOD(INT(AMASS+0.25),2).EQ.0) THEN
          IF (MOD(INT(ZNUM+0.25),2).EQ.0) THEN
            STAG = (EEXC - DENLO) / (DENHI - DENLO)
            IF (STAG.LE.0.0) STAG = 0.0
            IF (STAG.GE.1.0) STAG = 1.0
            IF (SPIN.LT.0.25) THEN
              FJ_LOC = FJ_LOC * (1.0 + 1.02 * (1.0 - STAG) )
            ELSEIF (MOD(INT(SPIN+0.25),2).EQ.0) THEN
              FJ_LOC = FJ_LOC * (1.0 + 0.227 * (1.0 - STAG) )
            ELSE
              FJ_LOC = FJ_LOC * (1.0 - 0.227 * (1.0 - STAG) )
            ENDIF
          ENDIF
        ENDIF
      ENDIF

C
C     "Parity-dependence" term
C
      IF (LDENP.EQ.0) THEN
        PARDEP=0.5
      ELSEIF (LDENP.EQ.2) THEN
        PAIRS=DENPEN(1)+DENPEN(2)/AMASS**DENPEN(3)    !see PRC67, 015803
      ELSEIF (LDENP.EQ.1) THEN                        !see PRC67, 015803
        IF (MOD(INT(AMASS+0.25),2).EQ.0) THEN
          IF (MOD(INT(ZNUM+0.25),2).EQ.0) THEN
            PAIRS= 1.34+75.22/AMASS**0.89             !E-E nucleus
          ELSE
            PAIRS=-0.90+75.22/AMASS**0.89             !O-O nucleus
          ENDIF
        ELSE
          IF (MOD(INT(ZNUM+0.25),2).EQ.0) THEN
            PAIRS=-0.08+75.22/AMASS**0.89             !E-O nucleus
          ELSE
            PAIRS=-0.42+75.22/AMASS**0.89             !O-E nucleus
          ENDIF
        ENDIF
      ENDIF
      IF (LDENP.GT.0) THEN
        IF (IPDEN.EQ.0) THEN 
          PARDEP=0.5*(1+1/(1+EXP(DENPC*(EEXC-PAIRS))))
        ELSE
          PARDEP=0.5*(1-1/(1+EXP(DENPC*(EEXC-PAIRS))))
        ENDIF
      ENDIF           
      IF (IPAR.EQ.1) PARDEP=1.0-PARDEP
c
      DENSITY=DENSITY*FJ_LOC*PARDEP
C
   31 CONTINUE
      RETURN
      END
c
C***********************************************************************
      FUNCTION ALD(EX,SPIN,IP)
C***********************************************************************
C
C     This subroutine provides cubic interpolation of values
C     of level density - based on the function AICC.
C     At the lowest excitations (very low density) the linear interpolation 
c     is used while at higher energies a cubic interpolation is adopted 
C                                                    Version from 16-MAY-09
C
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      DIMENSION XX(4),YY(4),A(4)
      REAL*4    YTMP
            INTEGER      NLMIN_TAB
      COMMON /LDTAB/ TABENLD(0:270),TABLD(0:270,0:MAXJC,0:1),
     *               TABTLD(0:270),NLD
      COMMON /LOGTAB/ TABENPSF_LOG(3,0:400),TABPSF_LOG(3,0:400),
     *                PSF_COEF(3,0:396,4),
     *                NLMIN_TAB(0:49,0:1)
C
      ALD = 0.0
      ISPIN = ISUBSC(SPIN)
c
c    Low level densities treated in a special way
c
      NLMIN = NLMIN_TAB(ISPIN,IP)
      IF (NLMIN.LT.0) GOTO 11
      IF (EX.LE.TABENLD(NLMIN)) GOTO 11   ! rho = 0
      IF (EX.LE.TABENLD(NLMIN+1)) THEN    ! linear interpolation for low rho
        EXPOM = (EX - TABENLD(NLMIN)) / 
     *          (TABENLD(NLMIN + 1) - TABENLD(NLMIN))
        ALD =TABLD(NLMIN+1,ISPIN,IP)-TABLD(NLMIN,ISPIN,IP)
        ALD = ALD * EXPOM
        GOTO 11
      ENDIF
      IF (EX.LE.TABENLD(NLMIN+2)) THEN    
       EXPOM = (EX - TABENLD(NLMIN + 1)) / 
     *         (TABENLD(NLMIN + 2) - TABENLD(NLMIN + 1))
       ALD=TABLD(NLMIN+2,ISPIN,IP)-TABLD(NLMIN+1,ISPIN,IP)
       ALD = TABLD(NLMIN+1,ISPIN,IP) + ALD * EXPOM
       GOTO 11
      ENDIF
c
      IF (EX.LE.TABENLD(2)) THEN
         K=0
         ELSE
         IF (EX.LE.TABENLD(NLD-1)) THEN
            ILO=3
            IHI=NLD-1
            DO WHILE (ILO.LT.IHI)
               IMID=(ILO+IHI)/2
               IF (EX.LE.TABENLD(IMID)) THEN
                  IHI=IMID
               ELSE
                  ILO=IMID+1
               ENDIF
            ENDDO
            K=ILO-3
            GO TO 1
            ELSE
            K=NLD-4
         ENDIF
      ENDIF
C
    1 exlog=log(ex)
      DO J=1,4
         if (TABENLD(K+J).lt.0.0) write(*,*)'A problem with NLD energy'
         XX(J)=log(TABENLD(K+J))
         YTMP = TABLD(K+J,ISPIN,IP)
         IF (YTMP.LE.0.0) YTMP = 1.E-2
         YY(J) = LOG(YTMP)
      ENDDO !J
      ALD=0.
      DO I=1,4
         A(I)=YY(I)
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)/(XX(J)-XX(I))
         ENDDO !J
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)*(XX(J)-EXLOG)
         ENDDO !J
         ALD=ALD+A(I)
      ENDDO !I
C
      ALD=exp(ALD)

   11 CONTINUE
      RETURN
      END
C
C***********************************************************************
      FUNCTION ATEMP(EX)
C***********************************************************************
C
C     Piecewise linear interpolation of nuclear temperature T(EX)
C     from Column 2 (TABTLD) of LDTAB.DAT when NOPTDE = 14.
C
      PARAMETER (MAXJC = 49)
      COMMON /LDTAB/ TABENLD(0:270),TABLD(0:270,0:MAXJC,0:1),
     *               TABTLD(0:270),NLD
      COMMON /BKTAB/  PSF_ENC(400,3),PSF_TTC(50,3),
     *                IBK_PSF(0:512,3),IBK_TT(0:64,3),IBK_LD(0:512),
     *                BK_E0(3),BK_SC(3),BK_T0(3),BK_TS(3),
     *                BK_L0,BK_LS,LBK_LD
      REAL*8          BK_E0,BK_SC,BK_T0,BK_TS,BK_L0,BK_LS
      INTEGER         IBK_PSF,IBK_TT,IBK_LD,NLT_BK
      LOGICAL         LBK_LD
C
      IF (EX.LE.0.0 .OR. NLD.LE.0) THEN
        ATEMP = 0.0
        RETURN
      ENDIF
      IF (EX.LE.TABENLD(1)) THEN
        IF (TABENLD(1).GT.0.0) THEN
          ATEMP = TABTLD(1) * (EX / TABENLD(1))
        ELSE
          ATEMP = TABTLD(1)
        ENDIF
        RETURN
      ENDIF
      IF (EX.GE.TABENLD(NLD)) THEN
        ATEMP = TABTLD(NLD)
        RETURN
      ENDIF
      IF (LBK_LD .AND. EX.EQ.EX) THEN
        ILO = NLT_BK(EX, TABENLD(1), NLD, 512, IBK_LD, BK_L0, BK_LS)
        GOTO 1601
      ENDIF
      ILO = 1
      IHI = NLD
      DO WHILE (ILO.LT.IHI-1)
        IMID = (ILO + IHI) / 2
        IF (EX.LE.TABENLD(IMID)) THEN
          IHI = IMID
        ELSE
          ILO = IMID
        ENDIF
      ENDDO
 1601 CONTINUE
      IF (TABENLD(ILO+1).GT.TABENLD(ILO)) THEN
        ATEMP = TABTLD(ILO) + (EX - TABENLD(ILO)) /
     *          (TABENLD(ILO+1) - TABENLD(ILO)) *
     *          (TABTLD(ILO+1) - TABTLD(ILO))
      ELSE
        ATEMP = TABTLD(ILO)
      ENDIF
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE OPEN_IT (INR,IRL)
C***********************************************************************
C
      COMMON /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS

C
      CHARACTER*3  EXTR,EXTS
      IERR=0
      IF (ISWBN.LT.0.OR.ISWBN.GT.2) THEN
         IERR=1
         RETURN
      ENDIF
      WRITE (EXTS,100) INR
      WRITE (EXTR,100) IRL
  100 FORMAT (I3.3)
      IF (ISWBN.EQ.1) THEN
         OPEN (UNIT=12,FILE='EVENTS.S'//EXTS//'.R'//EXTR,
     *        FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
         ELSE
         OPEN (UNIT=12,FILE='EVENTS.S'//EXTS//'.R'//EXTR,
     *        STATUS='UNKNOWN')
      ENDIF
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE DO_IT (NSTEP)
C***********************************************************************
C
      COMMON /monit/ ELQQ(0:126),SPQQ(0:126),IPQQ(0:126),ICQQ(126),
     &               DMQQ(126),WIQQ(0:126)
C$OMP THREADPRIVATE(/monit/)
      COMMON /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS

C
      IF (ISWBN.EQ.0) THEN
c         WRITE (12,*)   ' '
         WRITE (12,100) SPQQ(0),IPQQ(0),NSTEP,ELQQ(0)
         IF (ISWEL.EQ.1) WRITE (12,200) (ELQQ(J),J=1,NSTEP)
         IF (ISWSP.EQ.1) WRITE (12,300) (SPQQ(J),J=1,NSTEP)
         IF (ISWPA.EQ.1) WRITE (12,400) (IPQQ(J),J=1,NSTEP)
         IF (ISWIC.EQ.1) WRITE (12,400) (ICQQ(J),J=1,NSTEP)
         IF (ISWMX.EQ.1) WRITE (12,500) (DMQQ(J),J=1,NSTEP)
         IF (ISWWI.EQ.1) WRITE (12,500) (WIQQ(J),J=1,NSTEP)
      ELSEIF (ISWBN.EQ.2) THEN
         WRITE (12,600) SPQQ(0),IPQQ(0),NSTEP,ELQQ(0),
     *                  (ELQQ(J),J=1,NSTEP)
      ELSE
         WRITE (12) SPQQ(0),IPQQ(0),NSTEP,ELQQ(0)
         IF (ISWEL.EQ.1) WRITE (12) (ELQQ(J),J=1,NSTEP)
         IF (ISWSP.EQ.1) WRITE (12) (SPQQ(J),J=1,NSTEP)
         IF (ISWPA.EQ.1) WRITE (12) (IPQQ(J),J=1,NSTEP)
         IF (ISWIC.EQ.1) WRITE (12) (ICQQ(J),J=1,NSTEP)
         IF (ISWMX.EQ.1) WRITE (12) (DMQQ(J),J=1,NSTEP)
         IF (ISWWI.EQ.1) WRITE (12) (WIQQ(J),J=1,NSTEP)
      ENDIF
  100 FORMAT (F5.1,I3,I5,F11.5)
  200 FORMAT (2X,8F11.5)
  300 FORMAT (2X,8F11.1)
  400 FORMAT (2X,8I11)
  500 FORMAT (2X,8E11.3)
  600 FORMAT (F5.1,I3,I5,127F11.5)
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE CLOSE_IT
C***********************************************************************
C
      CLOSE (UNIT=12,STATUS='KEEP')
      RETURN
      END
C
c This file contain subroutines
c   LABELS - convert numeric value of give model to string label
c   WRITE_PARAMS - write down (to the output file) parameters concerning
c                  giant resonances, SP strength, Fermi liq. param., and
c                  values Ek0 and Eg0 (according to Kopecky)
c
c************************************************************************
      SUBROUTINE LABELS(NOPTDE,NOPTE1,NOPTM1,NOPTE2,tdens,tsfe1,tsfm1,
     &                      tsfe2)
c************************************************************************
        CHARACTER*(*) tdens,tsfe1,tsfm1,tsfe2
c
c    Conversion number of model to string label of the model
c
      tdens='???'
      tsfe1='???'
      tsfm1='???'
      tsfe2='???'

      IF (noptde.EQ.0) THEN
        tdens='CTF'
      ELSEIF (noptde.EQ.1) THEN
        tdens='BSFG'
      ELSEIF (noptde.EQ.14) THEN
        tdens='Tabulated NLD with tabulated temperature' //
     *        ' (Goriely RIPL)'
      ENDIF

      IF (nopte1.EQ.0) THEN
        tsfe1='SP'
      ELSEIF (nopte1.EQ.1) THEN
        tsfe1='BA'
      ELSEIF (nopte1.EQ.2) THEN
        tsfe1='TD-BA'
      ELSEIF (nopte1.EQ.3) THEN
        tsfe1='KMF+Chr'
      ELSEIF (nopte1.EQ.4) THEN
        tsfe1='KMF'
      ELSEIF (nopte1.EQ.5) THEN
        tsfe1='Chr'
      ELSEIF (nopte1.EQ.6) THEN
        tsfe1='Chr-phD'
      ELSEIF (nopte1.EQ.7) THEN
        tsfe1='SMLO'
      ELSEIF (nopte1.EQ.8) THEN
        tsfe1='8'
      ELSEIF (nopte1.EQ.9) THEN
        tsfe1='9'
      ELSEIF (nopte1.EQ.10) THEN
        tsfe1='10'
      ELSEIF (nopte1.EQ.11) THEN
        tsfe1='11'
      ENDIF

      IF (noptm1.EQ.0) THEN
        tsfm1='SP'
      ELSEIF (noptm1.EQ.1) THEN
        tsfm1='BA'
      ELSEIF (noptm1.EQ.2) THEN
        tsfm1='BAonSP'
      ELSEIF (noptm1.EQ.3) THEN
        tsfm1='Pow'
      ELSEIF (noptm1.EQ.4) THEN
        tsfm1='4'
      ENDIF

      IF (nopte2.EQ.0) THEN
        tsfe2='SP'
      ELSEIF (nopte2.EQ.1) THEN
        tsfe2='BA'
      ENDIF

      END

c************************************************************************
        SUBROUTINE WRITE_PARAMS(NAME)
c************************************************************************
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      CHARACTER*12   NAME
      CHARACTER*1    CHPAR
      COMMON /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     &               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     &               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     &               PAR_E1(3),DIPSLP,DIPZER
     &       /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
      COMMON /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
     &               NOPTCS,NLINc,CAPFR(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /vars/  factnrm
     &       /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
     &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
c
c       Writing down some inoput parameters to output file
c
      OPEN(9,FILE=NAME,STATUS='UNKNOWN')
      WRITE(9,*) 'Parameters of nuclear realizations:'
      WRITE(9,*) 'GDER parameters: Erez[MeV]   Width[MeV]   Sigma[mb]'
      DO i=1,NGIGE
        WRITE(9,111) ER(i),W0(i),SIG(i)
        ENDDO
  111   FORMAT('               ',2F11.2,F13.1)
      WRITE(9,*) 'GDMR parameters: Erez[MeV]   Width[MeV]   Sigma[mb]'
      DO i=1,NGIGM
        WRITE(9,111) ERM(i),WM0(i),SIGM(i)
      ENDDO
      WRITE(9,*) 'GQER parameters: Erez[MeV]   Width[MeV]   Sigma[mb]'
      DO i=1,NGIGE2
       WRITE(9,111) ERE(i),WE0(i),SIGE(i)
      ENDDO
      WRITE(9,112) DEG,DMG,QEL
  112 FORMAT(' SP strength:   E1:',E9.2,'   M1:',E9.2,'   E2:',E9.2)
      WRITE(9,114) TEMPER,EZERO,AMASS
  114 FORMAT(' Density parameters:  CTF:   Temper:',F5.3,'     E0:',
     &F6.3,'    Amass:',F5.0)
      WRITE(9,115) ASHELL,EONE,PAIRING
  115 FORMAT('                     BSFG:   Ashell:',F7.3,'  E1:',
     &F6.3,'  Pairing:',F6.3)
      WRITE(9,113) fermc,ek0,egzero,par_m1(1),par_m1(2),par_m1(3)
  113 FORMAT(' Fermi liq. par.:',F3.1,'            Ek0:',F4.1,'     Eg0:
     &',F4.1,'    PAR_M1(1):',E9.2,'  PAR_M1(2):',E9.2,
     &'  PAR_M1(3):',E9.2)
      IF (IPINc.EQ.0) THEN
       CHPAR='+'
      ELSE
       CHPAR='-'
      ENDIF
      WRITE(9,117) BN,SPINc(1),CHPAR,CAPFR(1)
  117 FORMAT(' Binding energy:',F7.4,'         Jpi:',F4.1,A1,'  Fract:',
     &F6.3)
      WRITE(9,116) factnrm
  116 FORMAT(' Normalization factor: ',E9.3)
      WRITE(9,*)'*******************************************************
     &***************'
      CLOSE(UNIT=9,STATUS='KEEP')
      END
C***********************************************************************
      SUBROUTINE LEVELS(SPC,IR)
C***********************************************************************
C
C    - no fluctuations at all LMODE=-1
C    - the Poisson distribution of neighbourhood level spacing is
C      assumed for LMODE=0 
C    - the Wigner distribution (with long-range correlations) is
C      assumed for LMODE=1 
c    - the "restricted Wigner distribution" - no long-range correlations
C
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
c
      integer   dekod,denum,delev,deparity
      INTEGER*4 IR,IRCON,IRCONc,NTOTAL
      REAL*4    EIN,EFI
C
      COMMON  /WID/  sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
      COMMON /MAINE/ eall,EIN,EFI,noptfl,BN,DELTA,NBIN,NTOTAL,iregi,
     *               noptcs,nlinc,capfr(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
     *               NEIGENVAL(0:MAXJC,0:1),LMODE
      COMMON  /GAU/   U,IFLAG
C$OMP THREADPRIVATE(/GAU/)
C
C     (Variables printed in lower case are not used in this subroutine)
C
      NTOTAL=0
      IFLAG=0

      IF (LMODE.EQ.0) THEN              !Poisson distribution
       DO IP=0,1
        SP=SPC-INT(SPC+.25)-1.
        DO J=0,MAXJC
         SP=SP+1.
         DO I=1,NBIN
          E=BN+DELTA/2.-FLOAT(I)*DELTA
          AVNL=DELTA*DENSITY(E,SP,IP)
          LEVCON(I,J,IP)=NPOISS(IR,AVNL)
          NTOTAL=NTOTAL+LEVCON(I,J,IP)
         ENDDO !I
        ENDDO !J
       ENDDO !IP
      ELSEIF (LMODE.EQ.1) THEN      !Wigner distribution with long-range correlations
       DO IP=0,1
        SP=SPC-FLOAT(INT(SPC+.25))-1.
        DO J=0,MAXJC
         SP=SP+1.
  	   X=10.*RAN0(IR)  ! Only to "randomize" the beginning of the sequence
	   K=1
    8	   IF (GOE_EIGEN_VAL(K,J,IP).LT.X) THEN
          K=K+1
          GOTO 8
	   ENDIF
         KLO=K
         DO I=NBIN,1,-1
          E=BN+DELTA/2.-FLOAT(I)*DELTA
          AVNL=DELTA*DENSITY(E,SP,IP)
	    X=X+AVNL
    9     IF (GOE_EIGEN_VAL(K,J,IP).LT.X) THEN 
           K=K+1
           GOTO 9
          ENDIF   
          KHI=K
          LEVCON(I,J,IP)=KHI-KLO        
          NTOTAL=NTOTAL+KHI-KLO
          KLO=KHI
         ENDDO !I
        ENDDO !J
       ENDDO !IP
      ELSEIF (LMODE.EQ.2) THEN              !Wigner distribution - no long-range correlations
       DO IP=0,1
        SP=SPC-INT(SPC+.25)-1.
         DO J=0,MAXJC
          DO I=0,NBIN
           LEVCON(I,J,IP)=0
          ENDDO  !I
          SP=SP+1.
          OMEGA=0.
          ALPHA=0.
          I=0
C
C         OMEGA is a random sample drawn from the Wigner distri-
C         bution; the expectation value of average distance between
C         neighbouring levels of a given spin and parity is assumed
C         to be equal to 1.
C
C           The constant 1.1283791671 is equal to "two divided by
C           square root of pi"
C
    1     RA=RAN0(IR)
          IF (RA.LE.0.) GO TO 1
          OMEGA=OMEGA+1.1283791671*SQRT(-ALOG(RA))
    3     IF (OMEGA.LT.ALPHA) THEN
            LEVCON(I,J,IP)=LEVCON(I,J,IP)+1
            NTOTAL=NTOTAL+1
            GO TO 1
          ELSE
            I=I+1
            IF (I.GT.NBIN) GO TO 2
            E=BN+DELTA/2.-FLOAT(I)*DELTA
            ALPHA=ALPHA+DENSITY(E,SP,IP)*DELTA
            GO TO 3
          ENDIF
    2     CONTINUE
         ENDDO !J
       ENDDO !IP
      ELSEIF (LMODE.EQ.-1) THEN             ! cumulative number of levels used; based on Wigner distribution but no randomnes
       DO IP=0,1
        SP=SPC-INT(SPC+.25)-1.
         DO J=0,MAXJC
          DO I=0,NBIN
           LEVCON(I,J,IP)=0
          ENDDO  !I
          SP=SP+1.
          OMEGA=0.
          ALPHA=0.
          I=0
C
   11     CONTINUE 
          OMEGA=OMEGA+1.0
   13     IF (OMEGA.LT.ALPHA) THEN
            LEVCON(I,J,IP)=LEVCON(I,J,IP)+1
            NTOTAL=NTOTAL+1
            GO TO 11
          ELSE
            I=I+1
            IF (I.GT.NBIN) GO TO 12
            E=BN+DELTA/2.-FLOAT(I)*DELTA
            ALPHA=ALPHA+DENSITY(E,SP,IP)*DELTA
            GO TO 13
          ENDIF
   12     CONTINUE
         ENDDO !J
       ENDDO !IP
      ELSEIF (LMODE.EQ.-2) THEN                  ! No fluctuations
       DO IP=0,1
        SP=SPC-INT(SPC+.25)-1.
        DO J=0,MAXJC
         SP=SP+1.
         DO I=1,NBIN
          E=BN+DELTA/2.-FLOAT(I)*DELTA
          AVNL=DELTA*DENSITY(E,SP,IP)
          LEVCON(I,J,IP)=NINT(AVNL)
          NTOTAL=NTOTAL+LEVCON(I,J,IP)
         ENDDO !I
        ENDDO !J
       ENDDO !IP


      ENDIF
C
C     NTOTAL is the total number of generated levels.
C
C     At this moment for each value of I the variable LEVCON(...) contains
C     the number of those levels of a particular spin and parity that fall
C     within the corresponding energy bin (whose width is DELTA).
C     The following DO-loops, however, convert this differential distribution
C     of level energies into a CUMULATIVE (i.e. integral) form. This simple
C     conversion leads to a significant increase of the speed of
C     function SEED.
C
    5 IAUX=0
      DO IP=0,1
       DO J=0,MAXJC
        LEVCON(0,J,IP)=IAUX
        DO I=1,NBIN
          LEVCON(I,J,IP)=LEVCON(I,J,IP)+LEVCON(I-1,J,IP)
        ENDDO 
        IAUX=LEVCON(NBIN,J,IP)
       ENDDO !J
      ENDDO !IP
      RETURN
      END
C
C***********************************************************************
      FUNCTION GOE_EIGEN_VAL(N,IS,IP)
C***********************************************************************
C
C     Eigenvalues were generated at the beginning 
C     by SUBROUTINE GENERATE_GOE_EIGEN_VAL(IR,N)
C
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      INTEGER*4 N
      COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
     *               NEIGENVAL(0:MAXJC,0:1),LMODE
C
      NN=NEIGENVAL(IS,IP)
      K=N-(N/NN)*NN
      GOE_EIGEN_VAL=EIGENVAL(K,IS,IP)+FLOAT((N/NN)*NN)
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE GENERATE_GOE_EIGEN_VAL(IR,N)
C***********************************************************************
C
C     Generates eigenvalues of random matrices - they are stored in
C     the COMMON /GOE/ and are used for generating level in the 
C     subroutine LEVELS via calling GOE_EIGEN_VAL
C
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      INTEGER   N,LDA
      PARAMETER (LDA=1001)
      REAL*8    A(LDA,N),RA(N),X,RAORD(N)
      INTEGER*4 IR
C
      COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
     *               NEIGENVAL(0:MAXJC,0:1),LMODE
      COMMON  /GAU/ U,IFLAG
C$OMP THREADPRIVATE(/GAU/)
C
      DATA      DPI  /3.141592653589793d0/ 
C
      IFLAG=0
      DO IP=0,1
       DO IS=0,MAXJC
        DO I=1,N
         A(I,I)=DBLE(GAUSS(IR)/SQRT(2.*FLOAT(N-1)))
         DO J=I+1,N
          A(I,J)=DBLE(GAUSS(IR)/SQRT(4.*FLOAT(N-1)))
          A(J,I)=A(I,J)
         ENDDO !J
        ENDDO !I
c        CALL DEVLSF(N,A,LDA,RA)
        CALL RSM1(LDA,N,A,RA)     
c        CALL SORTQQ(LOC(RA),N,SRT$REAL8)
        CALL SORT(N,RA,RAORD)
        DO I = 1,N
          RA(I) = RAORD(I)
        ENDDO

        M=-1
        DO L=1,N
         X=RA(L)	      
         IF (DABS(X).lt.0.9d0) THEN         !Only eigenval. in the middle are treated
          M=M+1
          EIGENVAL(M,IS,IP)=FLOAT(N-1)*
     *            (SNGL((X*DSQRT(1.d0-X**2)+DASIN(X))/DPI+.5d0))	    	      
         ENDIF		 			
        ENDDO !L
        NEIGENVAL(IS,IP)=M
        EIG0=EIGENVAL(0,IS,IP)
        DO MM=M,0,-1
         EIGENVAL(MM,IS,IP)=EIGENVAL(MM,IS,IP)-EIG0
        ENDDO !MM
       ENDDO !IS
      ENDDO !IP
      RETURN
      END

C
C*******************************************************************************
C
      SUBROUTINE SORT (N,X,XN)
C
C        ALGORITHM AS 304.8 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Sorts the N values stored in array X in ascending order
C
      INTEGER N
      DOUBLE PRECISION X(N),XN(N),XTEMP(N)
C
      INTEGER I, J, INCR
      DOUBLE PRECISION TEMP
C
	DO I = 1, N
       XTEMP(I) = X(I) 
	ENDDO
C
      INCR = 1
C
C        Loop : calculate the increment
C
   10 INCR = 3 * INCR + 1
      IF (INCR .LE. N) GOTO 10

C
C        Loop : Shell-Metzner sort
C
   20 INCR = INCR / 3
      I = INCR + 1
   30 IF (I .GT. N) GOTO 60
      TEMP = X(I)
      J = I
   40 IF (X(J - INCR) .LT. TEMP) GOTO 50
      X(J) = X(J - INCR)
      J = J - INCR
      IF (J .GT. INCR) GOTO 40
   50 X(J) = TEMP
      I = I + 1
      GOTO 30
   60 IF (INCR .GT. 1) GOTO 20
C
	DO I = 1, N
        XN(I) = X(I)
        X(I) = XTEMP(I) 
	ENDDO
C
      RETURN
      END
C
C*******************************************************************************
C
      SUBROUTINE rsm1(nm,n,a,w)
c
      integer n,nm,ierr
      integer k1,k2
c      double precision a(nm,n),w(n),fwork(1)
      double precision a(nm,n),w(n),fwork2(n),fwork1(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find all of the eigenvalues and some of the eigenvectors
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        m  the eigenvectors corresponding to the first m eigenvalues
c           are to be computed.
c           if m = 0 then no eigenvectors are computed.
c           if m = n then all of the eigenvectors are computed.
c
c     on output
c
c        w  contains all n eigenvalues in ascending order.
c
c        z  contains the orthonormal eigenvectors associated with
c           the first m eigenvalues.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat,
c           imtqlv and tinvit.  the normal completion code is zero.
c
c        fwork  is a temporary storage array of dimension 8*n.
c
c        iwork  is an integer temporary storage array of dimension n.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 10 * n
      if (n .gt. nm) go to 50
      k1 = 1
      k2 = k1 + n

c     .......... find eigenvalues only ..........

c      call  tred1(nm,n,a,w,fwork(k1),fwork(k2))
c      call  tqlrat(n,w,fwork(k2),ierr)

      call  tred1(nm,n,a,w,fwork1,fwork2)
      call  tqlrat(n,w,fwork2,ierr)

   50 return
      end
C
C*******************************************************************************
C
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      DOUBLE PRECISION D(N),E2(N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
C
C     This subroutine is a translation of the Algol procedure tqlrat,
C     Algorithm 464, Comm. ACM 16, 689(1973) by Reinsch.
C
C     This subroutine finds the eigenvalues of a symmetric
C     tridiagonal matrix by the rational QL method.
C
C     On input
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E2 contains the squares of the subdiagonal elements of the
C          input matrix in its last N-1 positions.  E2(1) is arbitrary.
C
C      On output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct and
C          ordered for indices 1,2,...IERR-1, but may not be
C          the smallest eigenvalues.
C
C        E2 has been destroyed.
C
C        IERR is set to
C          zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG for  DSQRT(A*A + B*B) .
C
C     Questions and comments should be directed to Burton S. Garbow,
C     Mathematics and Computer Science Div, Argonne National Laboratory
C
C     This version dated August 1987.
C     Modified by C. Moler to fix underflow/overflow difficulties,
C     especially on the VAX and other machines where epslon(1.0d0)**2
C     nearly underflows.  See the loop involving statement 102 and
C     the two statements just before statement 200.
C
C     ------------------------------------------------------------------
C

      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0D0
      T = 0.0D0
      E2(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
         H = DABS(D(L)) + DSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
         if (c .ne. 0.0d0) go to 105
C        Spliting tolerance underflowed.  Look for larger value.
         do 102 i = l, n
            h = dabs(d(i)) + dsqrt(e2(i))
            if (h .gt. t) t = h
  102    continue
         b = epslon(t)
         c = b * b
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = PYTHAG(P,1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
C           Avoid division by zero on next pass
            if (g .eq. 0.0d0) g = epslon(d(i))
            h = g * (p / r)
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
C
C*******************************************************************************
C
      subroutine tred1(nm,n,a,d,e,e2)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
c
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end
C
C*******************************************************************************
C
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
C
C*******************************************************************************
C
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
C

C***********************************************************************
      SUBROUTINE GERMS(IR)
C
C     This subroutine generates a random-generator random seed for
C     each individual level. The seeds obtained are stored in IRCON(K).
C     For a level of interest, characterized by {IP,SP,IB,IL}, the
C     corresponding seed can be fetched in a simple way. This is evident
C     from the body of the function SEED that performs such an operation.
C
C***********************************************************************
c
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      integer   dekod,denum,delev,deparity
      INTEGER*4 IR,IRCON,IRCONc,I,IG,K,N,NTOTAL,IAUX(8000000),
     &          NOLD,NDIF,NMGEN
      SAVE      IAUX
      REAL*4  EIN,EFI
      LOGICAL LPAIRED_GERMS
      COMMON /PGFLAG/ LPAIRED_GERMS
      INTEGER*8 SM_STATE, R64, C_GAMMA, SPLITMIX64
      PARAMETER (C_GAMMA = -7046029254386353131_8)
      INTEGER*4 ILEV, NTOT_ALL
C
      COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
      COMMON  /MAINE/ eall,EIN,EFI,noptfl,bn,delta,NBIN,NTOTAL,iregi,
     *                NOPTCS,NLINc,capfr(2)
C$OMP THREADPRIVATE(/MAINE/)

C
C     Default mode: SplitMix64 deterministic level seeding
C     (Used when --paired-germs CLI flag is not provided)
C
      IF (.NOT. LPAIRED_GERMS) THEN
        IRCONc(1) = IR
        SM_STATE = INT(IR, 8)
        IF (SM_STATE .LE. 0_8) SM_STATE = 123456789_8
        NTOT_ALL = NTOTAL + nddd
        DO 701 ILEV = 1, NTOT_ALL
          SM_STATE = SM_STATE + C_GAMMA
          R64 = SPLITMIX64(SM_STATE)
          IRCON(ILEV) = INT(MOD(IAND(R64, 9223372036854775807_8),
     &                          2147483646_8)) + 1
  701   CONTINUE
        IF (NOPTCS .NE. 1) THEN
          SM_STATE = SM_STATE + C_GAMMA
          R64 = SPLITMIX64(SM_STATE)
          IRCONc(2) = INT(MOD(IAND(R64, 9223372036854775807_8),
     &                         2147483646_8)) + 1
        ENDIF
        RETURN
      ENDIF
C
C     Legacy paired mode (when --paired-germs CLI flag is passed):
C     Seed for 1st capturing state
C
      k=1
      NOLD = 0
      IRCONc(k)=IR
C
C     IRCON(.) is seeded in a random fashion by seeds IR produced
C     consecutively by repeated calls of RAN(IR), see the body of the
C     routine. If two or more seeds are falling to the same word of
C     IRCON(.), then only the last of them is kept there.
C
      N = NTOTAL
      DO 1 K=1,NTOTAL+nddd
    1 IRCON(K)=0

      IF (N.GT.30000000) THEN
        KLAST = 1
        NMAXR = 0
    9   KLAST = KLAST + NMAXR 
        NMAXR = INT( N / 2) - 2
        DO K = KLAST, KLAST + NMAXR 
          DUMMY = RAN0(IR)
          IRCON(K) = IR
        ENDDO
        N = N - NMAXR
        IF (N.GT.15000000) GOTO 9 
      ENDIF

    4 CONTINUE
      IF (N.GT.15000000) THEN
        NRNDN = 100000000
      ELSE
        NRNDN = 10000000
      ENDIF
      DO 2 IG=1,NRNDN
      K=int(float(NTOTAL+nddd)*RAN0(IR))+1
      DUMMY=RAN0(IR)
    2 IRCON(K)=IR
C
C     Determine the number of not seeded sites N
C
      N=0
      DO 3 K=1,NTOTAL+nddd
      IF (IRCON(K).EQ.0) N=N+1
    3 CONTINUE

      NDIF = ABS(NOLD - N)
      NOLD = N
      IF ((N.GT.2000).AND.(NDIF.GT.1000)) GOTO 4
      IF (N.GT.2000) THEN
        IF (N.GT.8000000) THEN
          WRITE (*,*)
     &      'The level density is too high.'
          WRITE (*,*) ' '
          WRITE (*,*) 'The required dimension in subr GERMS is ', N
          STOP 'ERROR'
        ENDIF
      ENDIF
C
C     If N is too high (>2000), the WHOLE field is to be randomly
C     seeded with another handful of 10**6 seeds.
C
    7 IF (N.EQ.0) THEN
       IF (NOPTCS.NE.1) THEN
        DUMMY=RAN0(IR)
        k=2
        IRCONc(k)=IR         ! 2nd capture state is now seeded
       ENDIF
       RETURN
      ENDIF
C     If every site of the useful part of IRCON(.) is seeded,
C     the RETURN follows. If not, ONLY EMPTY WORDS  are additionally
C     seeded. A smaller number (100000) of seeds is used. Seeds
C     are again being distributed randomly. This new seeding to
C     a restricted area is to be repeated up to the point when ALL
C     words are seeded. After that the RETURN follows.
C
      N=0
      DO 5 K=1,NTOTAL+nddd
      IF (IRCON(K).NE.0) GOTO 5
      N=N+1
      IAUX(N)=K
    5 CONTINUE

      IF (N.LT.10000) THEN
        NMGEN = 100000
      ELSE
        NMGEN = 1000000
      ENDIF

      DO 6 IG=1,NMGEN
      I=int(float(N)*RAN0(IR))+1
      DUMMY=RAN0(IR)
    6 IRCON(IAUX(I))=IR
      GOTO 7
      END
C
C
C***********************************************************************
      SUBROUTINE WIDTHS (MODE,IPIN,SPIN,IBIN,ILIN,
     1             TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS)
C***********************************************************************
c    Diferrent treatment of primary transitions to low-lying levels
C
      INTEGER*4    MAXIS,MCACHE
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000,
     *              MCACHE = 4096            )   
      LOGICAL    renergy,SG_DONE(MAXBIN,4),CT_DONE(MAXBIN,0:1,3),
     *           LPAIRED_GERMS,LNO_WCACHE
      COMMON /PGFLAG/ LPAIRED_GERMS
      COMMON /CWOPT/  LNO_WCACHE
      REAL*4     SEFF1,SEFF2,GSQ1,GSQ2,GSQ,G1,G2,DENOM,CFAC1,CFAC2,
     *           GAUSS,GAUSS_SQ_ONE
      integer    dekod,denum,delev,deparity
      INTEGER*4  IRCON,IRCONc,ISCON,ISEED,SEEDS,ISDIS,ntotal,
     *           ILEV_W,ISLOT_W,C_TAG,C_IBIN,C_ISCON,C_ISDIS
      REAL       Im2Res
      REAL*4     EIN,EFI,SG_BIN(MAXBIN,4),CT_BIN(MAXBIN,0:1,3),
     *           C_EIN,C_TOTDIS,C_GADIS,C_STDIS
      REAL*8     TOTCON,STCON,GACON,DSUM_G1,DSUM_G2,
     *           C_TOTCON,C_GACON,C_STCON
      DIMENSION  TOTDIS(0:2),TOTCON(0:2),STCON(0:2,0:MAXBIN,-2:2,0:1),
     *           ISCON(0:2,0:MAXBIN,-2:2,0:1),GACON(0:2,-2:2,0:1),
     *           STDIS(0:2,0:30,-2:2,0:1),GADIS(0:2,-2:2,0:1),
     *           ISDIS(0:2,30,0:MAXJC,0:1),S(2),GG(2)
      COMMON  /WID/  sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd,dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
      COMMON /GAU/   u,IFLAG
C$OMP THREADPRIVATE(/GAU/)
      COMMON /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,ntotal,iregi,
     &               noptcs,nlinc,capfr(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /SFCE/  SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
     &               GE1DIS(2),GM1DIS(2),GE2DIS(2)
C$OMP THREADPRIVATE(/SFCE/)
      COMMON /enh/   sumg,ibmin,ibmax
     &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
     &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS
     &       /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
     &               CONVk(0:1,5,50)
      INTEGER*4 C_IFLAG,TAG_LOC
      REAL*4 C_U
      INTEGER*8 C_LOCK
      COMMON /WCACHE/ C_STCON(0:MAXBIN,-2:2,0:1,MCACHE),
     *                C_GACON(-2:2,0:1,MCACHE),C_TOTCON(MCACHE),
     *                C_ISCON(0:MAXBIN,-2:2,0:1,MCACHE),
     *                C_STDIS(0:30,-2:2,0:1,MCACHE),
     *                C_ISDIS(30,-2:2,0:1,MCACHE),
     *                C_GADIS(-2:2,0:1,MCACHE),C_TOTDIS(MCACHE),
     *                C_EIN(MCACHE),C_TAG(MCACHE),C_IBIN(MCACHE),
     *                C_IFLAG(MCACHE),C_U(MCACHE),
     *                C_LOCK(MCACHE)

C
C     Variables printed in lower case are not EXPLICITELY used in
C     this subroutine
C
C     If you use MODE=1 or 2 you have to set IBIN=1 and ILIN=1 !
C
      Q=0.0
      GE1=0.
      GM1=0.
      GE2=0.
      ILEV_W=0
      ISLOT_W=0

      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      IF (.NOT. LNO_WCACHE .AND. MODE.EQ.0 .AND. IREGI.EQ.0) THEN
        ILEV_W = LEVCON(IBIN-1,ISP,IPIN) + ILIN
        ISLOT_W = IAND(IEOR(ILEV_W,ISHFT(ILEV_W,-11)),MCACHE-1) + 1
        TAG_LOC = 0
C$OMP ATOMIC READ ACQUIRE
        TAG_LOC = C_TAG(ISLOT_W)
        IF (TAG_LOC.EQ.ILEV_W) THEN
C$OMP FLUSH
          EIN = C_EIN(ISLOT_W)
          TOTCON(0) = C_TOTCON(ISLOT_W)
          TOTDIS(0) = C_TOTDIS(ISLOT_W)
          IFLAG = C_IFLAG(ISLOT_W)
          U = C_U(ISLOT_W)
          DO IPFI=0,1
           DO ISPFI=ISP-2,ISP+2
            SPFI = SP + FLOAT(ISPFI)
            ISP_FI = INT(SPFI + .25)
            ISBS = NINT(SPFI + .25) - NINT(SPIN + .25)
            GACON(0,ISBS,IPFI) = C_GACON(ISBS,IPFI,ISLOT_W)
            GADIS(0,ISBS,IPFI) = C_GADIS(ISBS,IPFI,ISLOT_W)
            STCON(0,NBIN,ISBS,IPFI) = C_STCON(NBIN,ISBS,IPFI,ISLOT_W)
            IT = ITYPE(SPIN,IPIN,SPFI,IPFI)
            IF (IT.NE.0) THEN
             DO I=IBIN,NBIN
              STCON(0,I,ISBS,IPFI) = C_STCON(I,ISBS,IPFI,ISLOT_W)
             ENDDO
             DO I=IBIN+1,NBIN
              ISCON(0,I,ISBS,IPFI) = C_ISCON(I,ISBS,IPFI,ISLOT_W)
             ENDDO
             DO I=0,NDIS(ISP_FI,IPFI)
              STDIS(0,I,ISBS,IPFI) = C_STDIS(I,ISBS,IPFI,ISLOT_W)
             ENDDO
             DO I=1,NDIS(ISP_FI,IPFI)
              ISDIS(0,I,ISP_FI,IPFI) = C_ISDIS(I,ISBS,IPFI,ISLOT_W)
             ENDDO
            ENDIF
           ENDDO
          ENDDO
          RETURN
        ENDIF
      ENDIF

      IF (IREGI.EQ.0) THEN
        IF (MODE.NE.0) THEN
          EIN=BN
          Q=0.0
        ELSE
          NLEV=LEVCON(IBIN,ISUBSC(SPIN),IPIN)- 
     *         LEVCON(IBIN-1,ISUBSC(SPIN),IPIN)
          Q=1.-FLOAT(2*ILIN-1)/FLOAT(2*NLEV)
          EIN=BN-(FLOAT(IBIN)-Q)*DELTA
        ENDIF
        ISEED=SEEDS(MODE,SPIN,IPIN,ILIN,IBIN)
        IF (ISEED.LE.0) ISEED = 123456789
      ELSE
        EIN=EFI
      ENDIF
      SPAC=1./DENSITY(EIN,SPIN,IPIN)
      IF (SPAC.LE.0.0) SPAC=0.00001
      TOTCON(MODE)=0.D+0
      TOTDIS(MODE)=0.
c
      DO IS=-2,2
       DO IP=0,1
        STCON(MODE,NBIN,IS,IP)=0.D+0
        IF (MODE.EQ.0) THEN
         DO I=0,20
          STDIS(MODE,I,IS,IP)=0.
         ENDDO
        ENDIF
        GACON(MODE,IS,IP)=0.D+0
        GADIS(MODE,IS,IP)=0.
       ENDDO
      ENDDO
C
      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      IF (iregi.ne.1) THEN
       IF (MODE.GT.0) THEN
        DO I=IBIN+1,NBIN
         EG=(FLOAT(I-IBIN-1)+Q+.5)*DELTA
         DO ITT=1,4
          IF (ITT.NE.2) THEN
           SG_BIN(I,ITT)=SGAMMA(EG,EIN,ITT)
           SG_DONE(I,ITT)=.TRUE.
          ENDIF
         ENDDO
         DO MAEL0=0,1
          DO MUL0=1,3
           CT_BIN(I,MAEL0,MUL0)=AICCX(EG,1,MAEL0,MUL0)
           CT_DONE(I,MAEL0,MUL0)=.TRUE.
          ENDDO
         ENDDO
        ENDDO
       ELSE
        DO I=IBIN+1,NBIN
         SG_DONE(I,1)=.FALSE.
         SG_DONE(I,3)=.FALSE.
         SG_DONE(I,4)=.FALSE.
         DO MAEL0=0,1
          CT_DONE(I,MAEL0,1)=.FALSE.
          CT_DONE(I,MAEL0,2)=.FALSE.
          CT_DONE(I,MAEL0,3)=.FALSE.
         ENDDO
        ENDDO
       ENDIF
      ENDIF
      DO 1 IPFI=0,1
      DO 1 ISPFI= ISP-2, ISP+2
      SPFI = SP + FLOAT(ISPFI)
      ISP_FI = ISUBSC(SPFI)
      ISBS=NINT(SPFI+.25)-NINT(SPIN+.25)
C
      IT=ITYPE(SPIN,IPIN,SPFI,IPFI)
      IF (IT.EQ.0) GOTO 1
      IF (IT.EQ.2) THEN
         IT1=3
         IT2=4
         ELSE
         IT1=IT
         IT2=IT
      ENDIF
      MUL0=NINT(ABS(SPFI-SPIN))
      IF (MUL0.EQ.0) MUL0=1
      MUL1=MUL0+1
      IAUX=(-1)**(IPFI+IPIN+MUL0+1)
      IF (IAUX.EQ.-1) THEN
         MAEL0=0
      ELSE
         MAEL0=1
      ENDIF
      MAEL1=1-MAEL0
      STCON(MODE,IBIN,ISBS,IPFI)=0.D+0
      IF (iregi.eq.1) then
       ISEED=SEEDS(mode,spin,ipin,ilin,ibin)
       GOTO 24
      ENDIF
C
      DO I=IBIN+1,NBIN
       NL=LEVCON(I,ISP_FI,IPFI)-LEVCON(I-1,ISP_FI,IPFI)
       IF ((NOPTFL.EQ.1).AND.(NL.EQ.0)) THEN
        ISCON(MODE,I,ISBS,IPFI)=ISEED
        IFLAG=0
        STCON(MODE,I,ISBS,IPFI)=STCON(MODE,I-1,ISBS,IPFI)
        GOTO 25
       ENDIF
       EG=(FLOAT(I-IBIN-1)+Q+.5)*DELTA    ! Q=0.5 ... fixed; trans. between mid-bins
       Z=0.
       DO ITT=IT1,IT2
        IF (.NOT.SG_DONE(I,ITT)) THEN
         SG_BIN(I,ITT)=SGAMMA(EG,EIN,ITT)
         SG_DONE(I,ITT)=.TRUE.
        ENDIF
        S(ITT-IT1+1)=SG_BIN(I,ITT)
       ENDDO

       ISCON(MODE,I,ISBS,IPFI)=ISEED
       IF (NOPTFL.NE.1) THEN  !No fluctuation
         ZZ = 0.0
         GG(2)=0.0
         DO ITT=IT1,IT2
          GG(ITT-IT1+1) = S(ITT-IT1+1)
          ZZ = ZZ + GG(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
          DMIX2 = GG(2) / GG(1)
         ELSE
          DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,EG-EG,spin,ipin,dmix2,1)
         Z = Z + ZZ * FLOAT(NL) * (1 + alpha)
       ELSE
        IFLAG=0
        IF (.NOT.CT_DONE(I,MAEL0,MUL0)) THEN
         CT_BIN(I,MAEL0,MUL0)=AICCX(EG,1,MAEL0,MUL0)
         CT_DONE(I,MAEL0,MUL0)=.TRUE.
        ENDIF
        CT0=CT_BIN(I,MAEL0,MUL0)
        IF (IT.EQ.2) THEN
         IF (.NOT.CT_DONE(I,MAEL1,MUL1)) THEN
          CT_BIN(I,MAEL1,MUL1)=AICCX(EG,1,MAEL1,MUL1)
          CT_DONE(I,MAEL1,MUL1)=.TRUE.
         ENDIF
         CT1=CT_BIN(I,MAEL1,MUL1)
        ELSE
         CT1=0.0
        ENDIF
        SEFF1 = S(1)*(1.0 + CT0)
        IF (IT.EQ.2) THEN
          SEFF2 = S(2)*(1.0 + CT1)
        ELSE
          SEFF2 = 0.0
        ENDIF
        IF (MODE.NE.0 .AND. IVER.EQ.2) THEN
         DENOM = Re2Res(mode) + Im2Res(mode)
         CFAC1 = Im2Res(mode) / DENOM
         CFAC2 = Re2Res(mode) / ((1.0 + CORRI(mode)**2) * DENOM)
         IF (IT.EQ.2) THEN
          DSUM_G1 = 0.D+0
          DSUM_G2 = 0.D+0
          DO IL = 1, NL
           G1 = GAUSS(ISEED)
           G2 = GAUSS(ISEED) + CORRI(MODE)*G1
           GSQ1 = CFAC2*G2*G2 + CFAC1*G1*G1
           G1 = GAUSS(ISEED)
           G2 = GAUSS(ISEED) + CORRI(MODE)*G1
           GSQ2 = CFAC2*G2*G2 + CFAC1*G1*G1
           DSUM_G1 = DSUM_G1 + DBLE(GSQ1)
           DSUM_G2 = DSUM_G2 + DBLE(GSQ2)
          ENDDO
          Z = SNGL(DSUM_G1*DBLE(SEFF1) + DSUM_G2*DBLE(SEFF2))
         ELSE
          DSUM_G1 = 0.D+0
          DO IL = 1, NL
           G1 = GAUSS(ISEED)
           G2 = GAUSS(ISEED) + CORRI(MODE)*G1
           GSQ1 = CFAC2*G2*G2 + CFAC1*G1*G1
           DSUM_G1 = DSUM_G1 + DBLE(GSQ1)
          ENDDO
          Z = SNGL(DSUM_G1*DBLE(SEFF1))
         ENDIF
        ELSEIF (LPAIRED_GERMS) THEN
         IFLAG = 0
         IF (IT.EQ.2) THEN
          DSUM_G1 = 0.D+0
          DSUM_G2 = 0.D+0
          DO IL = 1, NL
           G1 = GAUSS(ISEED)
           G2 = GAUSS(ISEED)
           GSQ1 = G1*G1
           GSQ2 = G2*G2
           GG(1) = GSQ1*S(1)
           GG(2) = GSQ2*S(2)
           DSUM_G1 = DSUM_G1 + DBLE(GSQ1)
           DSUM_G2 = DSUM_G2 + DBLE(GSQ2)
          ENDDO
          Z = SNGL(DSUM_G1*DBLE(SEFF1) + DSUM_G2*DBLE(SEFF2))
         ELSE
          DSUM_G1 = 0.D+0
          DO IL = 1, NL
           G = GAUSS(ISEED)
           GSQ1 = G*G
           GG(1) = GSQ1*S(1)
           DSUM_G1 = DSUM_G1 + DBLE(GSQ1)
          ENDDO
          Z = SNGL(DSUM_G1*DBLE(SEFF1))
         ENDIF
        ELSE
         IF (IT.EQ.2) THEN
          DSUM_G1 = 0.D+0
          DSUM_G2 = 0.D+0
          CALL CHI2_SUM2(ISEED, NL, DSUM_G1, DSUM_G2)
          Z = SNGL(DSUM_G1*DBLE(SEFF1) + DSUM_G2*DBLE(SEFF2))
         ELSE
          DSUM_G1 = 0.D+0
          CALL CHI2_SUM1(ISEED, NL, DSUM_G1)
          Z = SNGL(DSUM_G1*DBLE(SEFF1))
         ENDIF
        ENDIF
       ENDIF ! NOPTFL
c
       STCON(MODE,I,ISBS,IPFI)=STCON(MODE,I-1,ISBS,IPFI)+
     *                         DBLE(Z*SPAC)
       IF (IT.EQ.1) THEN
         GE1=GE1+Z
       ELSEIF (IT.EQ.3) THEN
         GM1=GM1+Z
       ELSEIF (IT.EQ.4) THEN
         GE2=GE2+Z
       ELSEIF (IT.EQ.2) THEN          ! This is not fully correct (strength fluctuate)
         GM1=GM1+Z*S(1)/(S(1)+S(2))   ! But it is acceptable approximation
         GE2=GE2+Z*S(2)/(S(1)+S(2))
       ENDIF
   25  CONTINUE
      ENDDO !I
C
C     For primary transitions it is assumed that MODE=0. In such
C     a case it is understood that the values of the ACTUAL
C     subscripted variable that replaces the FICTIVE variable STDIS
C     are simply derived from input data without the need of Monte
C     Carlo simulation.
C
      IF (IPRIM.EQ.1) THEN
        IF (MODE.NE.0) GOTO 22    !!!!!!!! Compare to WIDTHS_R()
      ENDIF

   24 CONTINUE
      DO I=1,NDIS(ISP_FI,IPFI)
       EG=EIN-ENDIS(I,ISP_FI,IPFI)
       ISDIS(MODE,I,ISP_FI,IPFI)=ISEED   ! Needed only for mixing ratio
       IFLAG=0
       Z=0.
       DO ITT=IT1,IT2
        S(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
       ENDDO
       IF (NOPTFL.EQ.1) THEN
        IF (MODE.EQ.0) THEN
          ZZ = 0.0
          GG(2)=0.0
          DO ITT=IT1,IT2
           G=GAUSS(ISEED)
           GG(ITT-IT1+1) = G * G * S(ITT-IT1+1)
           ZZ = ZZ + GG(ITT-IT1+1)
          ENDDO
          IF (GG(1).GT.0.0) THEN
           DMIX2 = GG(2) / GG(1)
          ELSE
           DMIX2 = 0.0
          ENDIF
          alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,dmix2,1)
          Z = Z + ZZ * (1 + alpha)
        ELSE                    ! Primary transitions (the same fluctuation)
         ZZ = 0.0
         GG(2)=0.0
         DO ITT=IT1,IT2
           IF (IVER.EQ.2) THEN
            G1=GAUSS(ISEED)
            G2=GAUSS(ISEED)+CORRI(MODE)*G1
            GSQ=(Re2Res(mode)*G2*G2/(1+CORRI(mode)**2)+
     *           Im2Res(mode)*G1*G1)/(Re2Res(mode)+Im2Res(mode))
            GG(ITT-IT1+1) = GSQ * S(ITT-IT1+1)
           ELSE
            G=GAUSS(ISEED)
            GG(ITT-IT1+1) = G * G * S(ITT-IT1+1)
           ENDIF
           ZZ = ZZ + GG(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
           DMIX2 = GG(2) / GG(1)
         ELSE
           DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,dmix2,1)
         Z = Z + ZZ * (1 + alpha)
        ENDIF ! MODE
       ELSE
         ZZ = 0.0
         GG(2)=0.0
         DO ITT=IT1,IT2
          GG(ITT-IT1+1) = S(ITT-IT1+1)
          ZZ = ZZ + GG(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
          DMIX2 = GG(2) / GG(1)
         ELSE
          DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,0.0,spin,ipin,dmix2,1)
         Z = Z + ZZ * (1 + alpha)
       ENDIF ! NOPTFL
c
       IF (IT.EQ.1) THEN
         GE1=GE1+Z
       ELSEIF (IT.EQ.3) THEN
         GM1=GM1+Z
       ELSEIF (IT.EQ.4) THEN
         GE2=GE2+Z
       ELSEIF (IT.EQ.2) THEN           ! This is not fully correct (strength fluctuate)
         GM1=GM1+Z*S(1)/(S(1)+S(2))    ! But it is acceptable approximation
         GE2=GE2+Z*S(2)/(S(1)+S(2))
       ENDIF
c
       STDIS(MODE,I,ISBS,IPFI)=STDIS(MODE,I-1,ISBS,IPFI)+Z*SPAC     
      ENDDO

   22 TOTCON(MODE)=TOTCON(MODE)+STCON(MODE,NBIN,ISBS,IPFI)
      GACON(MODE,ISBS,IPFI)=TOTCON(MODE)
      TOTDIS(MODE)=TOTDIS(MODE)+
     &             STDIS(MODE,NDIS(ISP_FI,IPFI),ISBS,IPFI)
      GADIS(MODE,ISBS,IPFI)=TOTDIS(MODE)
C
    1 CONTINUE
      IF (.NOT. LNO_WCACHE) THEN
      IF (MODE.EQ.0 .AND. IREGI.EQ.0) THEN
        TAG_LOC = 0
C$OMP ATOMIC READ ACQUIRE
        TAG_LOC = C_TAG(ISLOT_W)
        IF (TAG_LOC.EQ.0) THEN
C$        CALL OMP_SET_LOCK(C_LOCK(ISLOT_W))
C$OMP ATOMIC READ ACQUIRE
          TAG_LOC = C_TAG(ISLOT_W)
          IF (TAG_LOC.EQ.0) THEN
            C_IBIN(ISLOT_W) = IBIN
            C_EIN(ISLOT_W) = EIN
            C_TOTCON(ISLOT_W) = TOTCON(0)
            C_TOTDIS(ISLOT_W) = TOTDIS(0)
            DO IPFI=0,1
             DO ISPFI=ISP-2,ISP+2
              SPFI = SP + FLOAT(ISPFI)
              ISP_FI = INT(SPFI + .25)
              ISBS = NINT(SPFI + .25) - NINT(SPIN + .25)
              C_GACON(ISBS,IPFI,ISLOT_W) = GACON(0,ISBS,IPFI)
              C_GADIS(ISBS,IPFI,ISLOT_W) = GADIS(0,ISBS,IPFI)
              C_STCON(NBIN,ISBS,IPFI,ISLOT_W) = STCON(0,NBIN,ISBS,IPFI)
              IT = ITYPE(SPIN,IPIN,SPFI,IPFI)
              IF (IT.NE.0) THEN
               DO I=IBIN,NBIN
                C_STCON(I,ISBS,IPFI,ISLOT_W) = STCON(0,I,ISBS,IPFI)
               ENDDO
               DO I=IBIN+1,NBIN
                C_ISCON(I,ISBS,IPFI,ISLOT_W) = ISCON(0,I,ISBS,IPFI)
               ENDDO
               DO I=0,NDIS(ISP_FI,IPFI)
                C_STDIS(I,ISBS,IPFI,ISLOT_W) = STDIS(0,I,ISBS,IPFI)
               ENDDO
               DO I=1,NDIS(ISP_FI,IPFI)
                C_ISDIS(I,ISBS,IPFI,ISLOT_W) = ISDIS(0,I,ISP_FI,IPFI)
               ENDDO
              ENDIF
             ENDDO
            ENDDO
            C_IFLAG(ISLOT_W) = IFLAG
            C_U(ISLOT_W) = U
C$OMP FLUSH
C$OMP ATOMIC WRITE RELEASE
            C_TAG(ISLOT_W) = ILEV_W
C$OMP FLUSH
          ENDIF
C$        CALL OMP_UNSET_LOCK(C_LOCK(ISLOT_W))
        ENDIF
      ENDIF
      ENDIF
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE WCACHE_INIT()
C***********************************************************************
      INTEGER*4 MAXBIN,MCACHE,C_TAG,C_IBIN,C_ISCON,C_ISDIS,I
      PARAMETER (MAXBIN = 700, MCACHE = 4096)
      REAL*4 C_EIN,C_TOTDIS,C_GADIS,C_STDIS
      REAL*8 C_TOTCON,C_GACON,C_STCON
      INTEGER*4 C_IFLAG
      REAL*4 C_U
      INTEGER*8 C_LOCK
      COMMON /WCACHE/ C_STCON(0:MAXBIN,-2:2,0:1,MCACHE),
     *                C_GACON(-2:2,0:1,MCACHE),C_TOTCON(MCACHE),
     *                C_ISCON(0:MAXBIN,-2:2,0:1,MCACHE),
     *                C_STDIS(0:30,-2:2,0:1,MCACHE),
     *                C_ISDIS(30,-2:2,0:1,MCACHE),
     *                C_GADIS(-2:2,0:1,MCACHE),C_TOTDIS(MCACHE),
     *                C_EIN(MCACHE),C_TAG(MCACHE),C_IBIN(MCACHE),
     *                C_IFLAG(MCACHE),C_U(MCACHE),
     *                C_LOCK(MCACHE)
      DO I=1,MCACHE
        C_TAG(I) = 0
        C_IBIN(I) = 999999
        C_IFLAG(I) = 0
        C_U(I) = 0.0
C$      CALL OMP_INIT_LOCK(C_LOCK(I))
      ENDDO
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE WCACHE_CLEAR()
C***********************************************************************
      INTEGER*4 MAXBIN,MCACHE,C_TAG,C_IBIN,C_ISCON,C_ISDIS,I
      PARAMETER (MAXBIN = 700, MCACHE = 4096)
      REAL*4 C_EIN,C_TOTDIS,C_GADIS,C_STDIS
      REAL*8 C_TOTCON,C_GACON,C_STCON
      INTEGER*4 C_IFLAG
      REAL*4 C_U
      INTEGER*8 C_LOCK
      COMMON /WCACHE/ C_STCON(0:MAXBIN,-2:2,0:1,MCACHE),
     *                C_GACON(-2:2,0:1,MCACHE),C_TOTCON(MCACHE),
     *                C_ISCON(0:MAXBIN,-2:2,0:1,MCACHE),
     *                C_STDIS(0:30,-2:2,0:1,MCACHE),
     *                C_ISDIS(30,-2:2,0:1,MCACHE),
     *                C_GADIS(-2:2,0:1,MCACHE),C_TOTDIS(MCACHE),
     *                C_EIN(MCACHE),C_TAG(MCACHE),C_IBIN(MCACHE),
     *                C_IFLAG(MCACHE),C_U(MCACHE),
     *                C_LOCK(MCACHE)
      DO I=1,MCACHE
        C_TAG(I) = 0
        C_IBIN(I) = 999999
        C_IFLAG(I) = 0
        C_U(I) = 0.0
      ENDDO
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE WCACHE_DESTROY()
C***********************************************************************
      INTEGER*4 MAXBIN,MCACHE,C_TAG,C_IBIN,I
      PARAMETER (MAXBIN = 700, MCACHE = 4096)
      INTEGER*4 C_IFLAG
      REAL*4 C_U,C_EIN,C_TOTDIS,C_GADIS,C_STDIS
      REAL*8 C_TOTCON,C_GACON,C_STCON
      INTEGER*8 C_LOCK
      COMMON /WCACHE/ C_STCON(0:MAXBIN,-2:2,0:1,MCACHE),
     *                C_GACON(-2:2,0:1,MCACHE),C_TOTCON(MCACHE),
     *                C_ISCON(0:MAXBIN,-2:2,0:1,MCACHE),
     *                C_STDIS(0:30,-2:2,0:1,MCACHE),
     *                C_ISDIS(30,-2:2,0:1,MCACHE),
     *                C_GADIS(-2:2,0:1,MCACHE),C_TOTDIS(MCACHE),
     *                C_EIN(MCACHE),C_TAG(MCACHE),C_IBIN(MCACHE),
     *                C_IFLAG(MCACHE),C_U(MCACHE),
     *                C_LOCK(MCACHE)
      DO I=1,MCACHE
C$      CALL OMP_DESTROY_LOCK(C_LOCK(I))
      ENDDO
      RETURN
      END
C
C
C
C***********************************************************************
      SUBROUTINE ONESTEP
     *   (MODE,IPIN,SPIN,IBIN,ILIN,TOTCON,STCON,GACON,ISCON,TOTDIS,
     *    STDIS,GADIS,ISDIS,IPFI,SPFI,IBFI,ILFI,DMIX2,SIGN,IR,IRX)
C***********************************************************************
C
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      logical    ponverze,ponverk,LPAIRED_GERMS
      COMMON /PGFLAG/ LPAIRED_GERMS
      integer    dekod,denum,delev,deparity,deaux,NLMIN_TAB
      REAL       Im2Res
      REAL*4     EIN,EFI
      REAL*8     AUX,AUX0,TOTCON,DRN,STCON,GACON,GAC,DZ,
     *           DSUM_G1,DSUM_G2
      INTEGER*4  ISEED,IR,IRX,ISDIS,ircon,irconc,ntotal,ISEEDEN,SEEDS
      REAL*4     SEFF1,SEFF2,STEP_Z,GSQ1,GSQ2,GSQ,G,G1,G2,Z
      REAL*4     DENOM,CFAC1,CFAC2,GAUSS,GAUSS_SQ_ONE
      DIMENSION  TOTDIS(0:2),TOTCON(0:2),STCON(0:2,0:MAXBIN,-2:2,0:1),
     *           ISCON(0:2,0:MAXBIN,-2:2,0:1),GACON(0:2,-2:2,0:1),
     *           STDIS(0:2,0:30,-2:2,0:1),GADIS(0:2,-2:2,0:1),
     *           ISDIS(0:2,30,0:MAXJC,0:1),S(2),GG(2)
C
      COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
      COMMON /GAU/   u,IFLAG
C$OMP THREADPRIVATE(/GAU/)
      COMMON /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,ntotal,iregi,
     *               noptcs,nlinc,capfr(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
     &               CONVk(0:1,5,50)
      COMMON /kon/   ponverze,ponverk
C$OMP THREADPRIVATE(/kon/)
      COMMON /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS
     &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
     &       /LOGTAB/ TABENPSF_LOG(3,0:400),TABPSF_LOG(3,0:400),
     *                PSF_COEF(3,0:396,4),
     *                NLMIN_TAB(0:49,0:1)
C
      IF (iregi.eq.2) GOTO 23
      SPAC=1./DENSITY(EIN,SPIN,IPIN)
      IF (SPAC.LE.0.0) SPAC=0.00001
C
    5 DRN=(DBLE(INT(DBLE(RAN0(IR))/1.D-4))+DBLE(RAN0(IR)))*1.D-4*
     *     (TOTCON(MODE)+DBLE(TOTDIS(MODE)))
      IF (DRN-TOTCON(MODE)) 1,2,2
C
C     Label #1 means that the transition ends in the continuum; this
C     can be learnt from IREGI=0.
C
    1 IREGI=0
      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      DO IPF=0,1
       DO ISPF = ISP-2, ISP+2
        SPF = SP + FLOAT(ISPF)
        ISBS=NINT(SPF+.25)-NINT(SPIN+.25)
        IT=ITYPE(SPIN,IPIN,SPF,IPF)
        IF (IT.NE.0) THEN
          GAC=GACON(MODE,ISBS,IPF)
          IF(DRN.LT.GAC) GOTO 4
        ENDIF
       ENDDO !SPF
      ENDDO !IPF
C
C     GOTO 5 statements are here for security reasons
C
      GOTO 5
C
C     Now IPFI,SPFI (& ISBS) are already determined:
C
    4 IPFI=IPF
      SPFI=SPF
      ISP_FI=ISUBSC(SPFI)
      DRN=DRN-GAC+STCON(MODE,NBIN,ISBS,IPFI)
      IF (IT.EQ.2) THEN
        IT1=3
        IT2=4
      ELSE
        IT1=IT
        IT2=IT
      ENDIF
C
C     Binary search for destination bin IBFI in cumulative STCON
C
      IF (IBIN+1.GT.NBIN) GOTO 5
      IF (STCON(MODE,NBIN,ISBS,IPFI).LE.DRN) GOTO 5
      ILO=IBIN+1
      IHI=NBIN
      DO WHILE (ILO.LT.IHI)
        IMID=(ILO+IHI)/2
        IF (STCON(MODE,IMID,ISBS,IPFI).GT.DRN) THEN
          IHI=IMID
        ELSE
          ILO=IMID+1
        ENDIF
      ENDDO
    7 IBFI=ILO
C
C     Now even IBFI is known. It remains to determine ILFI.
C
      ISEED=ISCON(MODE,IBFI,ISBS,IPFI)
      IFLAG=0
C
C     The appropriate ISEED is now fetched and IFLAG reset.
C     The algorithm for finding ILFI can start.
C
      EG=EIN-BN+(FLOAT(IBFI)-0.5)*DELTA
c      EG=(FLOAT(IBFI-IBIN))*DELTA	       !Energies between mid-bins
      NL=LEVCON(IBFI,ISP_FI,IPFI)-
     *   LEVCON(IBFI-1,ISP_FI,IPFI)
      AUX0=STCON(MODE,IBFI-1,ISBS,IPFI)
      DZ=0.D+0
      DO ITT=IT1,IT2
        S(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
      ENDDO
      IF (NL.GT.0) THEN
        CALL ALPH_PREP(EG,spfi,ipfi,0.0,spin,ipin,1,CT0,CT1)
      ELSE
        CT0 = 0.0
        CT1 = 0.0
      ENDIF
      SEFF1 = S(1)*(1.0 + CT0)
      IF (IT.EQ.2) THEN
        SEFF2 = S(2)*(1.0 + CT1)
      ELSE
        SEFF2 = 0.0
      ENDIF
      IF (NOPTFL.NE.1) THEN
        GG(1) = S(1)
        IF (IT.EQ.2) THEN
          GG(2) = S(2)
          STEP_Z = SEFF1 + SEFF2
        ELSE
          GG(2) = 0.0
          STEP_Z = SEFF1
        ENDIF
        DO IL=1,NL
          Z = FLOAT(IL) * STEP_Z
          AUX = AUX0 + DBLE(Z*SPAC)
          IF (AUX.GT.DRN) GOTO 9
        ENDDO
      ELSEIF (MODE.NE.0 .AND. IVER.EQ.2) THEN
        DENOM = Re2Res(mode) + Im2Res(mode)
        CFAC1 = Im2Res(mode) / DENOM
        CFAC2 = Re2Res(mode) / ((1.0 + CORRI(mode)**2) * DENOM)
        IF (IT.EQ.2) THEN
          DSUM_G1 = 0.D+0
          DSUM_G2 = 0.D+0
          DO IL=1,NL
            G1 = GAUSS(ISEED)
            G2 = GAUSS(ISEED) + CORRI(MODE)*G1
            GSQ1 = CFAC2*G2*G2 + CFAC1*G1*G1
            GG(1) = GSQ1*S(1)
            G1 = GAUSS(ISEED)
            G2 = GAUSS(ISEED) + CORRI(MODE)*G1
            GSQ2 = CFAC2*G2*G2 + CFAC1*G1*G1
            GG(2) = GSQ2*S(2)
            DSUM_G1 = DSUM_G1 + DBLE(GSQ1)
            DSUM_G2 = DSUM_G2 + DBLE(GSQ2)
            Z = SNGL(DSUM_G1*DBLE(SEFF1) + DSUM_G2*DBLE(SEFF2))
            AUX = AUX0 + DBLE(Z*SPAC)
            IF (AUX.GT.DRN) GOTO 9
          ENDDO
        ELSE
          DSUM_G1 = 0.D+0
          DO IL=1,NL
            G1 = GAUSS(ISEED)
            G2 = GAUSS(ISEED) + CORRI(MODE)*G1
            GSQ = CFAC2*G2*G2 + CFAC1*G1*G1
            GG(1) = GSQ*S(1)
            DSUM_G1 = DSUM_G1 + DBLE(GSQ)
            Z = SNGL(DSUM_G1*DBLE(SEFF1))
            AUX = AUX0 + DBLE(Z*SPAC)
            IF (AUX.GT.DRN) GOTO 9
          ENDDO
        ENDIF
      ELSEIF (LPAIRED_GERMS) THEN
        IFLAG = 0
        IF (IT.EQ.2) THEN
          DSUM_G1 = 0.D+0
          DSUM_G2 = 0.D+0
          DO IL=1,NL
            G1 = GAUSS(ISEED)
            G2 = GAUSS(ISEED)
            GSQ1 = G1*G1
            GSQ2 = G2*G2
            GG(1) = GSQ1*S(1)
            GG(2) = GSQ2*S(2)
            DSUM_G1 = DSUM_G1 + DBLE(GSQ1)
            DSUM_G2 = DSUM_G2 + DBLE(GSQ2)
            Z = SNGL(DSUM_G1*DBLE(SEFF1) + DSUM_G2*DBLE(SEFF2))
            AUX = AUX0 + DBLE(Z*SPAC)
            IF (AUX.GT.DRN) GOTO 9
          ENDDO
        ELSE
          DSUM_G1 = 0.D+0
          DO IL=1,NL
            G = GAUSS(ISEED)
            GSQ1 = G*G
            GG(1) = GSQ1*S(1)
            DSUM_G1 = DSUM_G1 + DBLE(GSQ1)
            Z = SNGL(DSUM_G1*DBLE(SEFF1))
            AUX = AUX0 + DBLE(Z*SPAC)
            IF (AUX.GT.DRN) GOTO 9
          ENDDO
        ENDIF
      ELSE
        IF (IT.EQ.2) THEN
          DSUM_G1 = 0.D+0
          DSUM_G2 = 0.D+0
          DO IL=1,NL
            CALL GAUSS_SQ_PAIR(ISEED, GSQ1, GSQ2)
            GG(1) = GSQ1*S(1)
            GG(2) = GSQ2*S(2)
            DSUM_G1 = DSUM_G1 + DBLE(GSQ1)
            DSUM_G2 = DSUM_G2 + DBLE(GSQ2)
            Z = SNGL(DSUM_G1*DBLE(SEFF1) + DSUM_G2*DBLE(SEFF2))
            AUX = AUX0 + DBLE(Z*SPAC)
            IF (AUX.GT.DRN) GOTO 9
          ENDDO
        ELSE
          DSUM_G1 = 0.D+0
          DO IL=1,NL
            GSQ = GAUSS_SQ_ONE(ISEED)
            GG(1) = GSQ*S(1)
            DSUM_G1 = DSUM_G1 + DBLE(GSQ)
            Z = SNGL(DSUM_G1*DBLE(SEFF1))
            AUX = AUX0 + DBLE(Z*SPAC)
            IF (AUX.GT.DRN) GOTO 9
          ENDDO
        ENDIF
      ENDIF
      GOTO 5
C
C     The resulting ILFI:
C
    9 ILFI=IL
C
      NLEV1=LEVCON(IBFI,ISP_FI,IPFI)-
     *      LEVCON(IBFI-1,ISP_FI,IPFI)
      Q1=1.-FLOAT(2*ILFI-1)/FLOAT(2*NLEV1)
c      ISEEDEN=SEEDS(0,SPFI,IPFI,ILFI,IBFI)
c      Q1=RAN0(ISEEDEN)
      EFI=BN-(FLOAT(IBFI)-Q1)*DELTA  ! Energy of the final level random in bin
      IF (efi.le.eall) IREGI=2
c
c       delta**2 and conversion - TODO - sign must be changed 
c
      IF (IT.EQ.2) THEN
        IF (GG(1).GT.0.0) THEN
          DMIX2 = GG(2) / GG(1)
        ELSE
          DMIX2 = 0.0
        ENDIF
        IFLAG=0
        GG(1) = GAUSS(ISEED)
        IF (GG(1).GT.0.) THEN
          SIGN=1.
        ELSE
          SIGN=-1.
        ENDIF
      ELSE
        DMIX2=0.
        SIGN=1.
      ENDIF
      dummy=RAN0(IRX)
      ALPHA=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,1)
      IF (dummy.GT.(ALPHA/(1+ALPHA))) THEN
        ponverze=.FALSE.
        ponverk=.FALSE.
      ELSE
        ALPHAK=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,2)
        ponverze=.TRUE.
        IF (dummy.LE.(ALPHAK/(1+ALPHAK))) THEN
          ponverk=.TRUE.
        ELSE
          ponverk=.FALSE.
        ENDIF
      ENDIF
      RETURN
C
    2 IREGI=1
      RN=SNGL(DRN-TOTCON(MODE))
      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      DO IPF=0,1
       DO ISPF = ISP-2, ISP+2
        SPF = SP + FLOAT(ISPF)
        ISBS=NINT(SPF+.25)-NINT(SPIN+.25)
        IT=ITYPE(SPIN,IPIN,SPF,IPF)
        IF (IT.NE.0) THEN
         IF (RN.LT.GADIS(MODE,ISBS,IPF)) GOTO 11
        ENDIF
       ENDDO !SPF
      ENDDO !IPF
      GOTO 5
   11 SPFI=SPF
      IPFI=IPF
      ISP_FI=ISUBSC(SPFI)
C
      RN=RN-GADIS(MODE,ISBS,IPFI)+
     *   STDIS(MODE,NDIS(ISP_FI,IPFI),ISBS,IPFI)
      DO I=1,NDIS(ISP_FI,IPFI)
        IF (RN.LT.STDIS(MODE,I,ISBS,IPFI)) GOTO 15
      ENDDO
      GOTO 5
   15 ILFI=I
      ibin=0                                        !!!!Why?
      EFI=ENDIS(ILFI,ISP_FI,IPFI)
      IF (EFI.LE.EALL) IREGI=2
c
c     delta**2 and conversion
c
      ISEED=ISDIS(MODE,ILFI,ISP_FI,IPFI)
      IFLAG=0
      IF (IT.EQ.2) THEN
         IT1=3
         IT2=4
         ELSE
         IT1=IT
         IT2=IT
      ENDIF
      EG=EIN-EFI
      DO ITT=IT1,IT2
        S(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
      ENDDO
      IF (NOPTFL.EQ.1) THEN
       IF (MODE.EQ.0) THEN
        DO ITT=IT1,IT2
         G=GAUSS(ISEED)
         GG(ITT-IT1+1)=G
        ENDDO
       ELSE                             !Primary transitions
        DO ITT=IT1,IT2
c         G1=GAUSS(ISEED)
c         G2=GAUSS(ISEED)+CORRI(MODE)*G1
c         GSQ=(Re2Res(mode)*G2*G2/(1+CORRI(mode)**2)+
c     *        Im2Res(mode)*G1*G1)/(Re2Res(mode)+Im2Res(mode))
c         GG(ITT-IT1+1)=sqrt(GSQ)
c         IF (G1.LT.0) GG(ITT-IT1+1)=-1.*GG(ITT-IT1+1)
         G=GAUSS(ISEED)
         GG(ITT-IT1+1)=G
        ENDDO
       ENDIF
      ENDIF
      IF (NOPTFL.NE.1) THEN
        GG(1)=1.
        GG(2)=1.
      ENDIF
      IF (IT.EQ.2) THEN
        IF (GG(1).NE.0. .AND. S(1).GT.0.0) THEN
          DMIX2 = GG(2)**2 * S(2) / (GG(1)**2 * S(1))
        ELSE
          DMIX2 = 1.E+6
        ENDIF
        IFLAG=0
        GG(1) = GAUSS(ISEED)
        IF (GG(1).GT.0.) THEN
          SIGN=1.
        ELSE
          SIGN=-1.
        ENDIF
      ELSE
        DMIX2=0.
        SIGN=1.
      ENDIF
      dummy=RAN0(IRX)
      ALPHA=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,1)
      IF (dummy.GT.(ALPHA/(1+ALPHA))) THEN
        ponverze=.FALSE.
        ponverk=.FALSE.
      ELSE
        ALPHAK=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,2)
        ponverze=.TRUE.
        IF (dummy.LE.(ALPHAK/(1+ALPHAK))) THEN
          ponverk=.TRUE.
        ELSE
          ponverk=.FALSE.
        ENDIF
      ENDIF
      RETURN
C
C
C     Note that this time ILFI means order number of an actual
C     level, not the order number of a level in an energy bin of
C     the energy continuum.
C
C     End in a part, where one knows every branching ratios
C
  23  EIN=EFI
      deaux=dekod(ilfi,isubsc(spfi),ipfi)
      rn=ran0(ir)*sall(deaux,denum(deaux))
      DO ia=1,denum(deaux)
        IF (rn.lt.sall(deaux,ia)) GOTO 24
      ENDDO
      GOTO 23
  24  ilfi=delev(deaux,ia)
      spfi=despin(deaux,ia)
      ipfi=deparity(deaux,ia)
      ibin=0
      efi=endis(ilfi,ISUBSC(spfi),ipfi)
      dmix2=deltx(deaux,ia)**2
      IF (deltx(deaux,ia).LT.0) THEN
        sign=-1.
      ELSE
        sign=1.
      ENDIF
c
C     conversion
c
      if (ponv(deaux,ia).gt.0.0) then
         dummy=ran0(ir)
         if (dummy.gt.ponv(deaux,ia)) then
           ponverze=.false.
           ponverk=.false.
         else
           ponverze=.true.
           if (dummy.le.ponvk(deaux,ia)) then
             ponverk=.true.
           else
             ponverk=.false.
           endif
         endif
      endif
      RETURN
      END
C
C***********************************************************************
      FUNCTION SEEDS(MODE,SP,IP,IL,IB)
C
C     This function yields the seed that is needed for the intrinsic
C     function RAN inside the subroutine WIDTHS immediately before
C     this subroutine starts to generate the values of the subscripted
C     variable STCON.
c
C***********************************************************************
C
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      INTEGER*4      SEEDS,IRCON,IRCONc
      integer        dekod,denum,delev,deparity
      COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
      COMMON /MAINE/  eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,ntotal,iregi,
     *               noptcs,nlinc,capfr(2)
C$OMP THREADPRIVATE(/MAINE/)
C
      IF (MODE.GT.0) THEN
        SEEDS=IRCONc(MODE)
      ELSE
        SEEDS=IRCON(LEVCON(IB-1,ISUBSC(SP),IP)+IL)
      ENDIF
C
      RETURN
      END
C
C*********************************************************************
      double precision function dran(ir)
c
c   this subroutine generates random values between 0.0 and 1.0 using
c   an integer seed
c   it is based on the imsl routine ggubs.
c   double precision version
c
      implicit double precision (a-h,o-z)
      integer*4 ir
      parameter(da=16807.d0,db=2147483647.d0,dc=2147483648.d0)
      ir=abs(mod(da*ir,db)+0.5d0)
      dran=dfloat(ir)/dc
      return
      end
c
c
      function sran(ir)
c
c   this subroutine generates random values between 0.0 and 1.0 using a 32-bit integer seed
c   it is based on the imsl routine ggubs.
c
      implicit double precision (a-h,o-z)
      integer*4 ir
      real*4 sran
      parameter(da=16807.d0,db=2147483647.d0,dc=2147483648.d0)
      ir=abs(mod(da*ir,db)+0.5d0)
      sran=float(ir)/dc
      return
      end
c
c
c   "Minimal" RND generator of Park and Miller - from Numerical Recipes
c
      function ran0(ir)
      integer*4 ir,IA,IM,IQ,IS,mask,K
      integer*8 P8
c      real*8    AM
      parameter (IA=16807,IM=2147483647,IQ=127773,IS=2836,
     *           mask=123459876)
c
c  other posible values of constants are:
c
c      parameter (IA=48271,IM=2147483647,IQ=44488,IS=3399,
c     *           mask=123459876)
c      parameter (IA=69621,IM=2147483647,IQ=30845,IS=23902,
c     *           mask=123459876)
c
c      ir=ieor(ir,mask)
   1  continue                   ! line added by MK
      AM=1./float(IM)
      IF (ir.GE.0) THEN
        P8=16807_8*INT(ir,8)
        P8=IAND(P8,2147483647_8)+ISHFT(P8,-31)
        IF (P8.GE.2147483647_8) P8=P8-2147483647_8
        ir=INT(P8,4)
      ELSE
        K=ir/IQ
        ir=IA*(ir-K*IQ)-IS*K
        if (ir.LT.0) ir=ir+IM
      ENDIF
c      ran0=AM*float(ir)
      ran0=AM*dble(ir)
c      ir=ieor(ir,mask)
       if (ran0.GE.1.0) goto 1   ! line added by MK
      return
      end
c
c  *********************************************************
c
c   "Minimal" RND generator of Park and Miller with Bays-Durham
c   shuffle and added safeguards - from Numerical Recipes
c
C     Dead code: function ran1 is not used anywhere in DICEBOX.
C     Commented out to remove static SAVE iv,iy variables.
Cc      function ran1(ir)
Cc      integer*4 ir,IA,IM,IQ,IS,NTAB,NDIV,K
Cc      parameter (IA=16807,IM=2147483647,IQ=127773,IS=2836,
Cc     *           NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
Cc
Cc  other posible values of constants are:
Cc
Cc      parameter (IA=48271,IM=2147483647,IQ=44488,IS=3399)
Cc      parameter (IA=69621,IM=2147483647,IQ=30845,IS=23902)
Cc
Cc      integer iv(NTAB)
Cc      save iv,iy
Cc      data iv /NTAB*0/, iy /0/
Cc
Cc      am=1./float(IM)
Cc      if (ir.LE.0.or.iy.EQ.0) then
Cc          ir=max(-ir,1)
Cc          do j=NTAB+8,1,-1
Cc             K=ir/IQ
Cc             ir=IA*(ir-K*IQ)-IS*K
Cc             if (ir.LT.0) ir=ir+IM
Cc             if (j.LE.NTAB) iv(j)=ir
Cc          enddo
Cc          iy=iv(1)
Cc      endif
Cc      K=ir/IQ
Cc      ir=IA*(ir-K*IQ)-IS*K
Cc      if (ir.LT.0) ir=ir+IM
Cc      j=1+iy/NDIV
Cc      iy=iv(j)
Cc      iv(j)=ir
Cc      ran1=min(AM*iy,RNMX)
Cc      return
Cc      end
c  ********************************************************
c
c   Congruental random number generator
c
        function ran2(ir)
        integer*4 ir,IA,IC,IM
      parameter (IA=7141,IM=259200,IC=54773)
c      parameter (IA=8121,IM=134456,IC=28411)

    1 continue
      ir=mod(ir*IA+IC,IM)
      if (ir.LT.0) goto 1
      ran2=float(ir)/float(IM)
      end

c***********************************************************************
      SUBROUTINE READ_EV(NAME)
c***********************************************************************

      INTEGER*4     MAXIS
      PARAMETER     (MAXJC  = 49,   !maximum spin generated in level scheme
     *               MAXBIN = 700,  !maximum # of bins in continuum
     *               MAXIS  = 76000000,
     *               MAXTPSF = 50     )   
      LOGICAL        LOGS,lopopgs,LEXIST_LD,LWARN_T
      CHARACTER*80   TITLE1,TITLE2,TITLE3
      CHARACTER*12   NAME,NAME_RES
      INTEGER*4      IR1,IR2,IR3,IR4,STEPS,IRCON,IRCONc,NTOTAL,IRINIT
      INTEGER        DEKOD,DENUM,DELEV,DEPARITY,NLMIN_TAB
      REAL       Im2Res,T_OLD,RECIP_T

      DIMENSION      ityp(2),qg(2),help(2),isbspin(2)      ! TNC
c
C     Hereafter, article 'CON' is related to a 'CONtinuum' part of
C     the level system, while 'DIS' -- to the 'DIScrete' part.
C     Suffixes 'p' and 's' are related to primary and secondary
C     transitions, respectively.
C
      COMMON /GAU/   U,IFLAG
C$OMP THREADPRIVATE(/GAU/)
      COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd,dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
     &       /incr/  nfilev,STEPS,enrgfin(18),buff(251),
     &               ncum,cumwidth(14),edo(18),eup(18)
      COMMON /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
     &               NOPTCS,NLINc,CAPFR(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     &               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     &               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     &               PAR_E1(3),DIPSLP,DIPZER
     &       /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
     &       /intr/  last,ir1,ir2,ir3,ir4,radwid,radwdi,
     &               fnoav,fnodi,rnn(18,14),
     &               rmult(18,126),rovfl(18),avea(18,14),erra(18,14),
     &               puvea(18,14),cova(18,18,14),aves(1:18,0:MAXBIN),
     &               errs(1:18,0:MAXBIN),
     &               diss(1:18,0:MAXBIN),poptlev(99),popslev(99)
     &       /readi/ nevents,NSREAL,NRL,ecrit,xrayk,xrayl,numlev,
     &               elowlev(99),elowsp(99),ilowip(99),
     &               STDISa(0:2,0:30,-2:2,0:1),RADW
     &       /grnd/  logs,kgs,LOpopGS,KpopGS
     &       /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
     &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
      COMMON /SFCE/  SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
     &               GE1DIS(2),GM1DIS(2),GE2DIS(2)
C$OMP THREADPRIVATE(/SFCE/)
      COMMON /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
     &               CONVk(0:1,5,50)
     &       /vars/  factnrm
     &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS
     &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
     &       /LOGTAB/ TABENPSF_LOG(3,0:400),TABPSF_LOG(3,0:400),
     *                PSF_COEF(3,0:396,4),
     *                NLMIN_TAB(0:49,0:1)
      COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
     *               NEIGENVAL(0:MAXJC,0:1),LMODE
     *       /denp/  LDENP,ZNUM,DENPC,DENPEN(3),IPDEN
      COMMON /LDTAB/ TABENLD(0:270),TABLD(0:270,0:MAXJC,0:1),
     *               TABTLD(0:270),NLD
      COMMON /BKTAB/  PSF_ENC(400,3),PSF_TTC(50,3),
     *                IBK_PSF(0:512,3),IBK_TT(0:64,3),IBK_LD(0:512),
     *                BK_E0(3),BK_SC(3),BK_T0(3),BK_TS(3),
     *                BK_L0,BK_LS,LBK_LD
      REAL*8          BK_E0,BK_SC,BK_T0,BK_TS,BK_L0,BK_LS
      INTEGER         IBK_PSF,IBK_TT,IBK_LD
      LOGICAL         LBK_LD
      COMMON /PSFTAB/TABENPSF(3,0:400),TABPSF(3,0:400),NPSF(3)
      COMMON /PSF2DTAB/ TABTPSF(3,MAXTPSF),TABPSF2D(3,0:400,MAXTPSF),
     *                  TABPSF2D_LOG(3,0:400,MAXTPSF),
     *                  PSF2D_COEF(3,0:396,4,MAXTPSF),
     *                  NTPSF(3),LWARN_T(3)
      COMMON /DISRAW/ ENRGF_RAW(99,20),SAL_RAW(99,20),ERRSAL_RAW(99,20),
     *                DESP_RAW(99,20),DLT_RAW(99,20),ALPHA_RAW(99,20),
     *                PRIM_RAW(99),ERRPRIM_RAW(99),IPAR_RAW(99,20)
     &       /rnd/   IRINIT(4,999)

      OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
C     OPEN (UNIT=4,FILE='AUXI.DAT')
C
C     User's alphanumeric titles:
C
      READ (5,100) TITLE1
      READ (5,100) TITLE2
      READ (5,100) TITLE3
  100 FORMAT (A80)
C
C     The regime of run:
C
      READ (5,*)
      READ (5,*) IVER,IPRIM
      READ (5,*)
      READ (5,*) ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,ISWIC,ISWMX,ISWWI,ISWLS
      READ (5,*)
      READ (5,*) NBIN, (IRINIT(I,1),I=1,4)
      IR1 = IRINIT(1,1) 
      IR2 = IRINIT(2,1)                   ! not really necessary
      IR3 = IRINIT(3,1)
      IR4 = IRINIT(4,1)
      READ (5,*)
      READ (5,*) NOPTFL,NOPTE1,NOPTM1,NOPTE2,NOPTDE,LMODE,LDENP
      READ (5,*)
      Read (5,*) NSREAL, NEVENTS, NRL
C
C     Giant Resonaces:
C
      READ (5,*)
      READ (5,*) NGIGE
      DO 15 I=1, NGIGE
        READ (5,*) ER(I),W0(I),SIG(I)
   15 CONTINUE
      READ (5,*)
      READ (5,*) NGIGM
      DO 16 I=1, NGIGM
        READ (5,*) ERM(I),WM0(I),SIGM(I)
   16 CONTINUE
      READ (5,*)
      READ (5,*) NGIGE2
      DO 17 I=1, NGIGE2
        READ (5,*) ERE(I),WE0(I),SIGE(I)
   17 CONTINUE
C
C     Other data needed for photon strength:
C
      READ (5,*)
      READ (5,*) DEG,DMG,QEL
      READ (5,*)
      READ (5,*) FERMC, EK0, EGZERO 
      READ (5,*)
      READ (5,*) (PAR_E1(I),I=1,3)
        DIPSLP=(1.0-PAR_E1(3))/(PAR_E1(2)-PAR_E1(1))
        DIPZER=1.0-DIPSLP*PAR_E1(2)
      READ (5,*)
      READ (5,*) (PAR_M1(I),I=1,3)
C
C     Data needed for level density formulas:
C
      READ (5,*)
      READ (5,*) ASHELL,EONE,TEMPER,EZERO,AMASS,ZNUM,PAIRING
      READ (5,*)
      READ (5,*) ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
      READ (5,*)
      READ (5,*) DENPC,(DENPEN(I),I=1,3)        !see PRC67, 015803 
      READ (5,*)
      READ (5,*) DENLO,DENHI,(DENPAR(I),I=1,4)
C
C     Data relating to the (thermal) neutron capturing state:
C
      READ (5,*)
      IF (IVER.EQ.0) THEN
        READ (5,*) BN,SPINCS,IPINC 
      ELSEIF (IVER.EQ.1) THEN
        READ (5,*) BN,SPINt,IPINC,CAPFR(1)
      ELSEIF (IVER.EQ.2) THEN
        READ (5,*) BN,SPINt,IPINC, En
        NAME_RES = 'RES_PAR.DAT'
        CALL READ_RES_PAR(NAME_RES,En,SPINt)
      ENDIF
      IF (IVER.EQ.0) THEN        ! IVER - resonances
        SPINc(1)=SPINCS
        NOPTCS=1
        NLINc=1
        CAPFR(1)=1.
      ELSEIF (IVER.EQ.1) THEN    ! TNC with fraction / for resonances already calculated in READ_RES_PAR()
       IF (SPINt.NE.0.) THEN
         NLINc=2
         NOPTCS=2
         spinc(1)=spint-0.5
         spinc(2)=spint+0.5
         capfr(2)=1.-capfr(1)
       ELSE
         NLINc=1
         NOPTCS=1
         spinc(1)=spint+0.5
         capfr(1)=1.
         capfr(2)=0.
       ENDIF
      ENDIF
C
C     Now read the tables of ICC-coefficients and other information
C     on electron conversion. The used notation:
C
C             XRAYK,XRAYL         -- electron binding energy for K-shell
C                                    (and L-shell) expressed in MeV
C             NENT,NENK           -- the number of electron energy points
C             NMU                 -- the highest multipolarity, e.g. 3 for
C                                    octupole (electric and magnetic)
C             CONVK(I,MU,K)       -- K-shell ICC for I-th type of radiation
C                                    (0 for electric, 1 for magnetic),
C                                    multipolarity MU and K-th electron
C                                    energy value
C             CONVT(I,M,K)        -- total ICC ... (as the CONVK)
C             ELENT(K),ELENK(K)   -- K-th value of electron energy
C
      read (5,*)
      read (5,*) xrayk,xrayl
      READ (5,*) NMU
      READ (5,*) NENT
      READ (5,*) (ELENT(K),K=1,NENT)
c      READ (5,*) (((CONVT(ITY,MU,K),MU=1,NMU),ITY=0,1),K=1,NENT)
      READ (5,*) (((CONVT(ITY,MU,K),ITY=0,1),MU=1,NMU),K=1,NENT)
      READ (5,*) NENK
      READ (5,*) (ELENK(K),K=1,NENK)
c      READ (5,*) (((CONVK(ITY,MU,K),MU=1,NMU),ITY=0,1),K=1,NENK)
      READ (5,*) (((CONVK(ITY,MU,K),ITY=0,1),MU=1,NMU),K=1,NENK)
      CALL AICC_PRECOMP()
C
C     Data related to the discrete levels (J, pi, Eexc, primary intensities
C     and branchings):
C
      do i=1,2
         GE1DIS(i)=0.
         GM1DIS(i)=0.
         GE2DIS(i)=0.
      enddo
      IPDEN = 0
c
      READ (5,*)
      READ (5,*) ECRIT
      EALL=ECRIT
      DELTA=(BN-ECRIT)/FLOAT(NBIN)
c      IF (IPRIM.GE.1) THEN   ! IPRIM
        READ (5,*)
        READ (5,*) FACTNRM
c      ENDIF

      auxi=0
      read (5,*)
      READ (5,*) numlev
      IF (NUMLEV.GT.99) STOP 'FATAL: NUMLEV > 99 exceeds DEKOD encoding'
      DO i=1,numlev
       IF (IPRIM.EQ.0) THEN
         READ (5,*) enrg,spfi,ipfi,denum(i)
         PRIM_RAW(i)=0.0
         ERRPRIM_RAW(i)=0.0
       ELSE 
         READ (5,*) enrg,spfi,ipfi,denum(i),prim,errprim
         PRIM_RAW(i)=prim
         ERRPRIM_RAW(i)=errprim
       ENDIF
       IF (denum(i).GT.20) THEN
         WRITE (*,*) 'number of de-exciting transitions is higher than',
     &               ' 20 (for level)',enrg
         STOP 'ERROR'  
       ENDIF
       IF (ndis(ISUBSC(spfi),ipfi).EQ.30) THEN
         WRITE (*,*) 'number of discrete levels with one spin/parity ',
     &               'is higher than 30',SPFI,IPFI
         STOP 'ERROR'  
       ENDIF
       ndis(ISUBSC(spfi),ipfi)=ndis(ISUBSC(spfi),ipfi)+1
       endis(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=enrg
       dekod(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=i
       IF (I.LE.5) THEN
         IPDEN = IPDEN + ipfi
       ENDIF
       IF (denum(i).GT.0) then
        DO k=1,denum(i)
         READ (5,*) enrgf,sal,errsal,desp,ipar,dlt,alpha
         ENRGF_RAW(i,k)=enrgf
         SAL_RAW(i,k)=sal
         ERRSAL_RAW(i,k)=errsal
         DESP_RAW(i,k)=desp
         IPAR_RAW(i,k)=ipar
         DLT_RAW(i,k)=dlt
         ALPHA_RAW(i,k)=alpha
         if (sal.LE.0.0) then
           WRITE(*,*)'intensity must be >0 (between levels)',enrg,enrgf
           STOP 'ERROR'
         endif 
         DO il=1,ndis(ISUBSC(desp),ipar)
           IF (ABS(endis(il,ISUBSC(desp),ipar)-enrgf).LE.1e-5) THEN
            delev(i,k)=il
            despin(i,k)=desp
            deparity(i,k)=ipar
            deltx(i,k)=dlt
            dmixsq=dlt**2
            if (alpha.LE.1e-6) alpha=ALPH_TOT(enrg,spfi,ipfi,enrgf,
     *                            desp,ipar,dmixsq,1)
            alphak=ALPH_TOT(enrg,spfi,ipfi,enrgf,desp,ipar,dmixsq,2)
            if (alphak.GT.alpha) alphak=alpha
            ponv(i,k)=alpha/(1+alpha)
            ponvk(i,k)=alphak/(1+alpha)
            sall(i,k)=sall(i,k-1)+sal*(1+alpha)
           ENDIF
         ENDDO
         if (delev(i,k).EQ.0) then
           WRITE(*,*) 'Wrong final level (Eini, Efin)',enrg,enrgf
           STOP 'ERROR'
         endif
        ENDDO
       ENDIF
       elowlev(i)=enrg
       IF(elowlev(i).EQ.0) THEN
        LOpopGS=.TRUE.
        KpopGS=i
       ENDIF
       elowsp(i)=spfi
       ilowip(i)=ipfi
       IF (IPDEN.GT.2) THEN 
         IPDEN = 1
       ELSE
         IPDEN = 0
       ENDIF
c
c     Intensities of primary transitions
c
       IF (IPRIM.GE.1) THEN
        DO ILINc=1,NLINc
          isbspin(ilinc)=NINT(spfi+.25)-NINT(spinc(ILINc)+.25)
          ityp(ilinc)=ITYPE(spinc(ilinc),ipinc,spfi,ipfi)
        ENDDO
        IF (NLINc.EQ.1) ityp(2)=0
        IF ((ityp(1).EQ.0).AND.(ityp(2).EQ.0)) THEN
          IF (prim.GT.0) THEN
            write(*,*) "No possible transition for primary feeding"
            write(*,*) "for level",i," with spin/parity",spfi,ipfi
            STOP
          ENDIF
        ELSEIF ((ityp(1).NE.0).AND.(ityp(2).EQ.0)) THEN
         STDISa(1,ndis(ISUBSC(spfi),ipfi),isbspin(1),ipfi)=
     &         STDISa(1,ndis(ISUBSC(spfi),ipfi)-1,isbspin(1),ipfi)
     &         +prim*factnrm/capfr(1)
         IF (ityp(1).EQ.1) ge1dis(1)=ge1dis(1)+prim*factnrm/capfr(1)
         IF ((ityp(1).EQ.3).OR.(ityp(1).EQ.2)) gm1dis(1)=gm1dis(1)
     &                                +prim*factnrm/capfr(1)
         IF (ityp(1).EQ.4) ge2dis(1)=ge2dis(1)+prim*factnrm/capfr(1)
        ELSEIF ((ityp(1).EQ.0).AND.(ityp(2).NE.0)) THEN
          STDISa(2,ndis(ISUBSC(spfi),ipfi),isbspin(2),ipfi)=
     &         STDISa(2,ndis(ISUBSC(spfi),ipfi)-1,isbspin(2),ipfi)
     &         +prim*factnrm/capfr(2)
         IF (ityp(2).EQ.1) ge1dis(2)=ge1dis(2)+prim*factnrm/capfr(2)
         IF ((ityp(2).EQ.3).OR.(ityp(2).EQ.2)) gm1dis(2)=gm1dis(2)
     &                                +prim*factnrm/capfr(2)
         IF (ityp(2).EQ.4) ge2dis(2)=ge2dis(2)+prim*factnrm/capfr(2)
        ELSEIF ((ityp(1).NE.0).AND.(ityp(2).NE.0)) THEN
         DO ILINc=1,NLINc
           IF (prim.GT.0) THEN
             QG(ilinc)=sgamma(bn-enrg,bn,ityp(ilinc))/
     *                 density(bn,spinc(ilinc),ipinc)/(prim*factnrm)
           ELSE
             QG(ilinc)=0.
           ENDIF
         ENDDO
         qgsum=0.
         DO ilinc=1,nlinc
           qgsum=qgsum+qg(ilinc)
         ENDDO
         DO ilinc=1,nlinc
            if (qgsum.GT.0) then
              help(ilinc)=QG(ilinc)/qgsum*prim*factnrm/capfr(ilinc)
            else
              help(ilinc)=0.
            endif
           ENDDO
         DO ilinc=1,nlinc
           STDISa(ilinc,ndis(ISUBSC(spfi),ipfi),isbspin(ilinc),ipfi)=
     &       STDISa(ilinc,ndis(ISUBSC(spfi),ipfi)-1,isbspin(ilinc),ipfi)
     &       +help(ilinc)
           IF (ityp(ilinc).EQ.1)
     &       ge1dis(ilinc)=ge1dis(ilinc)+help(ilinc)
           IF ((ityp(ilinc).EQ.3).OR.(ityp(ilinc).EQ.2))
     &       gm1dis(ilinc)=gm1dis(ilinc)+help(ilinc)
           IF (ityp(ilinc).EQ.4)
     &       ge2dis(ilinc)=ge2dis(ilinc)+help(ilinc)
         ENDDO
        ENDIF        ! ityp
c      
       ENDIF ! IPRIM
c
      ENDDO  ! i - NUMLEV
c
      DO 31 ipfi=0,1
       sp=spinc(1)-int(spinc(1)+.25)
       DO 31 ispfi = 0, 8
        spfi = sp + float(ispfi)
        levdis(ISUBSC(spfi),ipfi)=auxi+ndis(ISUBSC(spfi),ipfi)
        auxi=levdis(ISUBSC(spfi),ipfi)
   31   CONTINUE
      nddd=auxi
c
c      KGS=0
c      DO I=1,NFILEV
c        READ (5,*) ENRGFIN(I)
c        IF (ENRGFIN(I).EQ.0.) THEN
c          LOGS=.TRUE.
c          KGS=I
c        ENDIF
c        EDO(I)=ENRGFIN(I)-1.E-6
c        EUP(I)=ENRGFIN(I)+1.E-6
c      ENDDO
c
      CLOSE (UNIT=5,STATUS='KEEP')
c
c     Tabulated level density (Goriely / Temperature-coupled LDTAB.DAT)
c
c
c     Tabulated level density (Goriely: NOPTDE=11, or Goriely with T column: NOPTDE=14)
c
      IF (NOPTDE.EQ.11 .OR. NOPTDE.EQ.14) THEN
        NAME = "LDTAB.DAT"
        OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
        READ(5,*) 
        READ(5,*) NLD, SPACRES,SPINi,IPINi,corAlpLD,corDelLD
        READ(5,*) 
        TABENLD(0) = 0.0
        TABTLD(0) = 0.0
        DO IEN = 1, NLD
          READ(5,*) TABENLD(IEN),TABTLD(IEN),DUMMY,DUMMY,DUMMY,
     *              (TABLD(IEN,J,0),J=0,MAXJC) 
        ENDDO
        READ(5,*) 
        DO IEN = 1, NLD
          READ(5,*) TABENLD(IEN),TABTLD(IEN),DUMMY,DUMMY,DUMMY,
     *              (TABLD(IEN,J,1),J=0,MAXJC) 
        ENDDO
        CLOSE(5)
c
c     "Normalization to experimental level spacing" - SPINt/IPINt not needed at the moment
c      PRC78,064307(2008): \rho(u,J,p)=\exp(alpha x \sqrt(U-delta)) * \rho(U-delta,J,p)
c
c  Normalization is not "correct": PRC78,064307(2008) gives \rho(u,J,p)=\exp(alpha x \sqrt(U-delta)) * \rho(U-delta,J,p)
c  parameters in RIPL correcton file are (likely) alpha = ctable and delta = ptable
c
        DO I = 1, NLD
          IF (NOPTDE.EQ.14 .AND. corAlpLD.NE.0.0 .AND.
     *        TABENLD(I).GT.0.0 .AND. TABTLD(I).GT.0.0) THEN
            T_OLD = TABTLD(I)
            RECIP_T = 1.0 / T_OLD + corAlpLD / (2.0 * SQRT(TABENLD(I)))
            IF (RECIP_T .GT. 0.0) THEN
              TABTLD(I) = 1.0 / RECIP_T
            ENDIF
          ENDIF
          FSPAC = EXP(corAlpLD * sqrt(TABENLD(I)))
          DO J = 0, MAXJC
           DO K = 0, 1
             TABLD(I,J,K) = TABLD(I,J,K) * FSPAC
           ENDDO
          ENDDO
          TABENLD(I)=TABENLD(I)+corDelLD
        ENDDO
      ENDIF 
c
c     Tabulated level density (Kawano: LDTAB_K.DAT) - Exact Prague 1.0.1
c
      IF (NOPTDE.EQ.12) THEN
        NAME = "LDTAB_K.DAT"
        OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
        READ(5,*) 
        READ(5,*) NLD
        READ(5,*) 
        TABENLD(0) = 0.0
        TABTLD(0) = 0.0
        DO IEN = NLD, 1, -1
          READ(5,*) TABENLD(IEN),DUMMY,
     *              (TABLD(IEN,J,0),J=0,9) 
        ENDDO
        CLOSE(5)
        DO IEN = 1, NLD
          DO J = 0, 9
            TABLD(IEN,J,0) = TABLD(IEN,J,0) / 2.0 ! Toshihiko gives total NLD (p-independent)
            TABLD(IEN,J,1) = TABLD(IEN,J,0) 
          ENDDO
          DO J = 10, MAXJC
            TABLD(IEN,J,0) = 0.0
            TABLD(IEN,J,1) = 0.0
          ENDDO
        ENDDO
      ENDIF
c
c     Precompute NLMIN_TAB (first positive entry of TABLD for each spin/parity)
c
      DO KPAR = 0, 1
        DO JSPIN = 0, MAXJC
          NLMIN_TAB(JSPIN,KPAR) = -1
        ENDDO
      ENDDO
      IF (NOPTDE.EQ.11 .OR. NOPTDE.EQ.12 .OR. NOPTDE.EQ.14) THEN
        DO KPAR = 0, 1
          DO JSPIN = 0, MAXJC
            DO I = 1, NLD
              IF (TABLD(I,JSPIN,KPAR).GT.0.0) THEN
                NLMIN_TAB(JSPIN,KPAR) = I - 1
                GOTO 209
              ENDIF
            ENDDO
 209        CONTINUE
          ENDDO
        ENDDO
        CALL BKT_BUILD(TABENLD(1), NLD, 512, IBK_LD, BK_L0, BK_LS,
     *                 LBK_LD)
      ENDIF
c
c     Tabulated PSF (1D and 2D temperature-dependent formats)
c
      DO I = 1, 3
        NTPSF(I) = 0
        LWARN_T(I) = .FALSE.
      ENDDO
      IF ((NOPTE1.EQ.11).OR.(NOPTE1.EQ.12)) THEN
        NAME = "PSFE1.DAT"
        CALL READ_PSF_FILE(1, NAME)
      ENDIF

      IF ((NOPTM1.GE.11).AND.(NOPTM1.LE.13)) THEN
        NAME = "PSFM1.DAT"
        CALL READ_PSF_FILE(2, NAME)
      ENDIF

      IF ((NOPTE2.EQ.11).OR.(NOPTE2.EQ.12)) THEN
        NAME = "PSFE2.DAT"
        CALL READ_PSF_FILE(3, NAME)
      ELSEIF ((NOPTE2.GE.13).AND.(NOPTE2.LE.15)) THEN
        WRITE(*,*) 'FATAL: Unsupported NOPTE2 = ', NOPTE2
        STOP 'FATAL: Unsupported NOPTE2'
      ENDIF
C
      IF ((NTPSF(1).GE.2).OR.(NTPSF(2).GE.2).OR.(NTPSF(3).GE.2)) THEN
        IF ((NOPTDE.EQ.11).OR.(NOPTDE.EQ.12)) THEN
          WRITE(*,*) 'WARNING: Tabular NLD (NOPTDE=', NOPTDE,
     *               ') lacks explicit temperature for 2D PSF.'
          WRITE(*,*) '         Falling back to Fermi-gas temperature ',
     *               'or T clamping.'
        ENDIF
      ENDIF
C
      WRITE(*,*) ' Inputs have been successfully loaded.'

      RETURN
      END
c
c***********************************************************************
      SUBROUTINE READ_INT(NAME)
c***********************************************************************
c
      INTEGER*4    MAXIS
      REAL       Im2Res
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      CHARACTER*12   NAME
      INTEGER*4      IR1,IR2,IR3,IR4,STEPS,IRCON,IRCONc,NTOTAL
      INTEGER        DEKOD,DENUM,DELEV,DEPARITY,NLMIN_TAB
c
C     Hereafter, article 'CON' is related to a 'CONtinuum' part of
C     the level system, while 'DIS' -- to the 'DIScrete' part.
C     Suffixes 'p' and 's' are related to primary and secondary
C     transitions, respectively.
C
      COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),irconc(2),deltx(99,20)
     &       /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     &               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     &               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     &               PAR_E1(3),DIPSLP,DIPZER
     &       /readi/ nevents,NSREAL,NRL,ecrit,xrayk,xrayl,numlev,
     &               elowlev(99),elowsp(99),ilowip(99),
     &               STDISa(0:2,0:30,-2:2,0:1),RADW
     &       /intr/  last,ir1,ir2,ir3,ir4,radwid,radwdi,
     &               fnoav,fnodi,rnn(18,14),
     &               rmult(18,126),rovfl(18),avea(18,14),erra(18,14),
     &               puvea(18,14),cova(18,18,14),aves(1:18,0:MAXBIN),
     &               errs(1:18,0:MAXBIN),
     &               diss(1:18,0:MAXBIN),poptlev(99),popslev(99)
     &       /incr/  nfilev,STEPS,enrgfin(18),buff(251),
     &               ncum,cumwidth(14),edo(18),eup(18)
      COMMON /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
     &               NOPTCS,NLINc,CAPFR(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
     &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
      COMMON /SFCE/  SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
     &               GE1DIS(2),GM1DIS(2),GE2DIS(2)
C$OMP THREADPRIVATE(/SFCE/)
      COMMON /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
     &               CONVk(0:1,5,50)
     &       /vars/  factnrm
      COMMON /GAU/   U,IFLAG
C$OMP THREADPRIVATE(/GAU/)
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
     &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS
     &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
     &       /LOGTAB/ TABENPSF_LOG(3,0:400),TABPSF_LOG(3,0:400),
     *                PSF_COEF(3,0:396,4),
     *                NLMIN_TAB(0:49,0:1)
      COMMON /DISRAW/ ENRGF_RAW(99,20),SAL_RAW(99,20),ERRSAL_RAW(99,20),
     *                DESP_RAW(99,20),DLT_RAW(99,20),ALPHA_RAW(99,20),
     *                PRIM_RAW(99),ERRPRIM_RAW(99),IPAR_RAW(99,20)
c
      iflag=0
      DO J=0,MAXJC
        DO K=0,1
          NDIS(J,K)=0
        ENDDO
      ENDDO
c
      DO i=1,numlev
       enrg=elowlev(i)
       spfi=elowsp(i)
       ipfi=ilowip(i)
       prim=PRIM_RAW(i)
       errprim=ERRPRIM_RAW(i)
       IF (ndis(ISUBSC(spfi),ipfi).GE.30)
     *   STOP 'FATAL: NDIS >= 30 in READ_INT'
       ndis(ISUBSC(spfi),ipfi)=ndis(ISUBSC(spfi),ipfi)+1
       endis(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=enrg
       dekod(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=i
       IF (denum(i).GT.0) then
        DO k=1,denum(i)
         enrgf=ENRGF_RAW(i,k)
         sal=SAL_RAW(i,k)
         errsal=ERRSAL_RAW(i,k)
         desp=DESP_RAW(i,k)
         ipar=IPAR_RAW(i,k)
         dlt=DLT_RAW(i,k)
         alpha=ALPHA_RAW(i,k)
         DO il=1,ndis(ISUBSC(desp),ipar)
           IF (ABS(endis(il,ISUBSC(desp),ipar)-enrgf).LE.1.0E-5) THEN
            delev(i,k)=il
            despin(i,k)=desp
            deparity(i,k)=ipar
            deltx(i,k)=dlt
            dmixsq=dlt**2
            if (alpha.LE.1e-6) alpha=ALPH_TOT(enrg,spfi,ipfi,enrgf,
     *                            desp,ipar,dmixsq,1)
            alphak=ALPH_TOT(enrg,spfi,ipfi,enrgf,desp,ipar,dmixsq,2)
            if (alphak.GT.alpha) alphak=alpha
            ponv(i,k)=alpha/(1+alpha)
            ponvk(i,k)=alphak/(1+alpha)
   51       sal=sal+errsal*GAUSS(IR4)             !!! not fixed
            if (sal.LE.0) goto 51
            sall(i,k)=sall(i,k-1)+sal*(1+alpha)
           ENDIF
         ENDDO
        ENDDO
       ENDIF
      ENDDO
c
      RETURN
      END
c
c***********************************************************************
      SUBROUTINE READ_AGAIN(NAME)
c***********************************************************************
c
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      CHARACTER*12   NAME
      INTEGER*4      IR1,IR2,IR3,IR4,STEPS,IRCON,IRCONc,NTOTAL
      INTEGER        DEKOD,DENUM,DELEV,DEPARITY,NLMIN_TAB
      REAL       Im2Res
      DIMENSION      QG(2),ITYP(2),ISBSPIN(2)

      CHARACTER*80  xxx2

c
C     Hereafter, article 'CON' is related to a 'CONtinuum' part of
C     the level system, while 'DIS' -- to the 'DIScrete' part.
C     Suffixes 'p' and 's' are related to primary and secondary
C     transitions, respectively.
C
      COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),irconc(2),deltx(99,20)
     &       /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     &               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     &               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     &               PAR_E1(3),DIPSLP,DIPZER
     &       /readi/ nevents,NSREAL,NRL,ecrit,xrayk,xrayl,numlev,
     &               elowlev(99),elowsp(99),ilowip(99),
     &               STDISa(0:2,0:30,-2:2,0:1),RADW
     &       /intr/  last,ir1,ir2,ir3,ir4,radwid,radwdi,
     &               fnoav,fnodi,rnn(18,14),
     &               rmult(18,126),rovfl(18),avea(18,14),erra(18,14),
     &               puvea(18,14),cova(18,18,14),aves(1:18,0:MAXBIN),
     &               errs(1:18,0:MAXBIN),
     &               diss(1:18,0:MAXBIN),poptlev(99),popslev(99)
     &       /incr/  nfilev,STEPS,enrgfin(18),buff(251),
     &               ncum,cumwidth(14),edo(18),eup(18)
      COMMON /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
     &               NOPTCS,NLINc,CAPFR(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
     &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
      COMMON /SFCE/  SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
     &               GE1DIS(2),GM1DIS(2),GE2DIS(2)
C$OMP THREADPRIVATE(/SFCE/)
      COMMON /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
     &               CONVk(0:1,5,50)
      COMMON /GAU/   U,IFLAG
C$OMP THREADPRIVATE(/GAU/)
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
     &       /vars/  factnrm
     &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS
     &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
     &       /LOGTAB/ TABENPSF_LOG(3,0:400),TABPSF_LOG(3,0:400),
     *                PSF_COEF(3,0:396,4),
     *                NLMIN_TAB(0:49,0:1)
      COMMON /DISRAW/ ENRGF_RAW(99,20),SAL_RAW(99,20),ERRSAL_RAW(99,20),
     *                DESP_RAW(99,20),DLT_RAW(99,20),ALPHA_RAW(99,20),
     *                PRIM_RAW(99),ERRPRIM_RAW(99),IPAR_RAW(99,20)
c
c         Preparing of some variables
c
      DO K=0,1
        DO L=0,20
          DO M=-2,2
            DO MM=0,2
              STDISa(MM,L,M,K)=0.
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO J=0,MAXJC
        DO K=0,1
          NDIS(J,K)=0
        ENDDO
      ENDDO
      DO i=1,2
        GE1DIS(i)=0.
        GM1DIS(i)=0.
        GE2DIS(i)=0.
      ENDDO
c
      DO i=1,numlev
       enrg=elowlev(i)
       spfi=elowsp(i)
       ipfi=ilowip(i)
       prim=PRIM_RAW(i)
       errprim=ERRPRIM_RAW(i)
       IF (ndis(ISUBSC(spfi),ipfi).GE.30)
     *   STOP 'FATAL: NDIS >= 30 in READ_AGAIN'
       ndis(ISUBSC(spfi),ipfi)=ndis(ISUBSC(spfi),ipfi)+1
       endis(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=enrg
       dekod(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=i
       DO ILINc=1,NLINc
         isbspin(ilinc)=NINT(spfi+.25)-NINT(spinc(ILINc)+.25)
         ityp(ilinc)=ITYPE(spinc(ilinc),ipinc,spfi,ipfi)
       ENDDO
       IF ((ityp(1).EQ.0).AND.(ityp(2).EQ.0)) THEN
         IF (prim.GT.0) THEN
           write(*,*) "No possible transition"
           STOP
         ENDIF
       ELSEIF ((ityp(1).NE.0).AND.(ityp(2).EQ.0)) THEN
        STDISa(1,ndis(ISUBSC(spfi),ipfi),isbspin(1),ipfi)=
     &        STDISa(1,ndis(ISUBSC(spfi),ipfi)-1,isbspin(1),ipfi)
     &        +prim*factnrm/capfr(1)
        IF (ityp(1).EQ.1) ge1dis(1)=ge1dis(1)+prim*factnrm/capfr(1)
        IF ((ityp(1).EQ.3).OR.(ityp(1).EQ.2)) gm1dis(1)=gm1dis(1)
     &                               +prim*factnrm/capfr(1)
        IF (ityp(1).EQ.4) ge2dis(1)=ge2dis(1)+prim*factnrm/capfr(1)
       ELSEIF ((ityp(1).EQ.0).AND.(ityp(2).NE.0)) THEN
         STDISa(2,ndis(ISUBSC(spfi),ipfi),isbspin(2),ipfi)=
     &        STDISa(2,ndis(ISUBSC(spfi),ipfi)-1,isbspin(2),ipfi)
     &        +prim*factnrm/capfr(2)
        IF (ityp(2).EQ.1) ge1dis(2)=ge1dis(2)+prim*factnrm/capfr(2)
        IF ((ityp(2).EQ.3).OR.(ityp(2).EQ.2)) gm1dis(2)=gm1dis(2)
     &                               +prim*factnrm/capfr(2)
        IF (ityp(2).EQ.4) ge2dis(2)=ge2dis(2)+prim*factnrm/capfr(2)
       ELSEIF ((ityp(1).NE.0).AND.(ityp(2).NE.0)) THEN
         QGnfl=0.
         DO ILINc=1,NLINc
          if (prim.GT.0.) then
            QGnfl=QGnfl+capfr(ilinc)*sgamma(bn-enrg,bn,ityp(ilinc))/
     *            density(bn,spinc(ilinc),ipinc)/(prim*factnrm)
          else
            QGnfl=1.
          endif
         ENDDO
  151    CONTINUE
         DO ILINc=1,NLINc
          qg(ilinc)=0.
          iflag=0
          IF (prim.GT.0) THEN
            IF (IVER.LE.1) THEN
              G1=GAUSS(ir4)
              QG(ilinc)=sgamma(bn-enrg,bn,ityp(ilinc))*G1*G1/
     *                  density(bn,spinc(ilinc),ipinc)/(prim*factnrm)
            ELSEIF (IVER.EQ.2) THEN
              G1=GAUSS(ir4)
              G2=GAUSS(ir4)+CORRI(ilinc)*G1
              GSQ=(Re2Res(ilinc)*G2*G2/(1+CORRI(ilinc)**2)+
     *            Im2Res(ilinc)*G1*G1)/(Re2Res(ilinc)+Im2Res(ilinc))
              QG(ilinc)=sgamma(bn-enrg,bn,ityp(ilinc))*GSQ/
     *                  density(bn,spinc(ilinc),ipinc)/(prim*factnrm)
            ENDIF
          ELSE
            QG(ilinc)=0.
          ENDIF
         ENDDO
         IF (((abs(qg(1)*capfr(1)+qg(2)*capfr(2)-QGnfl))/QGnfl.LE.0.01)
     &      .OR.(prim.EQ.0)) THEN
          DO ilinc=1,nlinc
           STDISa(ilinc,ndis(ISUBSC(spfi),ipfi),isbspin(ilinc),ipfi)=
     &       STDISa(ilinc,ndis(ISUBSC(spfi),ipfi)-1,isbspin(ilinc),ipfi)
     &       +prim*factnrm*QG(ilinc)/QGnfl
           IF (ityp(ilinc).EQ.1) ge1dis(ilinc)=ge1dis(ilinc)+
     &                                  prim*factnrm*QG(ilinc)/QGnfl
           IF ((ityp(ilinc).EQ.3).OR.(ityp(ilinc).EQ.2))
     &          gm1dis(ilinc)=gm1dis(ilinc)+prim*factnrm*QG(ilinc)/QGnfl
           IF (ityp(ilinc).EQ.4) ge2dis(ilinc)=ge2dis(ilinc)+
     &                                  prim*factnrm*QG(ilinc)/QGnfl
          ENDDO
          GOTO 152
         ENDIF
         GOTO 151
       ENDIF
  152  CONTINUE
      ENDDO
      RETURN
      END
c
c***********************************************************************
      SUBROUTINE READ_RES_PAR(NAME,En,SPINt)
c***********************************************************************
c
      REAL           Im2Res
      CHARACTER*12   NAME
      DIMENSION      ERES(50,2),GRED(50,2),GGAM(50,2),NRES(2),
     &               gf(2),sg(2),rho(2),RI(2)
     
      COMMON /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
     &               NOPTCS,NLINc,CAPFR(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
     &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
     &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
C
C     Data about resonances
C
      OPEN(UNIT=15,FILE=NAME,STATUS='OLD')
      READ(15,*)
      READ(15,*) NRES(1)
      DO l=1,NRES(1)
       READ(15,*) ERES(l,1),GRED(l,1),GGAM(l,1)
      ENDDO
      IF (SPINt.NE.0.) THEN
        READ(15,*)
        READ(15,*) NRES(2)
        DO l=1,NRES(2)
         READ(15,*) ERES(l,2),GRED(l,2),GGAM(l,2)
        ENDDO
      ENDIF
      CLOSE(15)
c
c     Auxiliary variables connected with resonances
c
      DO j=1,2
       Re2Res(j)=0.
       Im2Res(j)=0.
       RI(j)=0.
       sg(j)=0.
      ENDDO
      IF (SPINt.NE.0) THEN
       jmax=2
       gf(1)=0.5*(2*(spint-0.5)+1)/(2*spint+1)
       gf(2)=0.5*(2*(spint+0.5)+1)/(2*spint+1)
      ELSE
       jmax=1
       gf(1)=0.5*(2*(spint+0.5)+1)/(2*spint+1)
       gf(2)=gf(1)
      ENDIF
c
      DO j=1,jmax
       DO l=1,NRES(j)
        gtot=(GGAM(l,j)+sqrt(En)/2/gf(j)*GRED(l,j))
        denom=((En-ERES(l,j))**2+1.e-6*gtot**2./4)
        rnoma=(ERES(l,j)-En)*sqrt(GRED(l,j)/2./gf(j))
        rnomb=1.e-3*gtot/2.*sqrt(GRED(l,j)/2./gf(j))
        Re2Res(j)=Re2Res(j)+(rnoma/denom)**2
        Im2Res(j)=Im2Res(j)+(rnomb/denom)**2
        RI(j)=RI(j)+rnoma/denom*rnomb/denom
        sigma=GGAM(l,j)*GRED(l,j)/2/((En-ERES(l,j))**2+
     *        1e-6/4*(GGAM(l,j)+sqrt(En)/2/gf(j)*GRED(l,j))**2)
        sg(j)=sg(j)+sigma
       ENDDO
       rho(j)=RI(j)/(sqrt(Re2Res(j))*sqrt(Im2Res(j)))
       CORRI(j)=rho(j)/sqrt(1-rho(j)**2)
      ENDDO
      capfr(1)=sg(1)/(sg(1)+sg(2))
      capfr(2)=1.-capfr(1)
c
c      write(*,*) 'Capturing spin fractions:',capfr(1),capfr(2)
c
      IF (jmax.EQ.1) THEN
       NLINc=1
       NOPTCS=1
       spinc(1)=spint+0.5
       spinc(2)=spint+1.5
      ELSEIF (capfr(2).LT.0.001) THEN    ! One per cent contribution is the minimum
       capfr(1)=1.                      ! for consideration of capture spin
       capfr(2)=0.
       NLINc=1
       NOPTCS=1
       spinc(1)=spint-0.5
       spinc(2)=spint+0.5
      ELSEIF (capfr(1).LT.0.001) THEN
       capfr(1)=1.
       capfr(2)=0.
       NLINc=1
       NOPTCS=1
       spinc(1)=spint+0.5
       spinc(2)=spint+1.5
       Re2Res(1)=Re2Res(2)
       Im2Res(1)=Im2Res(2)
       CORRI(1)=CORRI(2)
      ELSE
       NLINc=2
       NOPTCS=2
       spinc(1)=spint-0.5
       spinc(2)=spint+0.5
      ENDIF

      RETURN
      END
c
c***********************************************************************
c                                             version 2.54
c
c   - changed meaning of some of NOPTE1 (3,12-17) in this version
c
c Functions TERM(Excit_En)
c           SGAM(Egam,Eini,Ityp)
c
C***********************************************************************
      FUNCTION SGAMMA(EGAM,EINI,ITYP)
c
c   Photon strengths for E1, M1+E2, M1 and E2 transitions. "Photon
c   strength" does not mean "photon strength function" here, but
c   photon strength function multiplied by EGAM**(2*L+1) and, in
c   the case of M1+E2 transitions, summed over both XL-components.
c
c     ITYP is equal to 1, 2, 3 or 4 (see ITYPE function)
c     EINI is the initial state energy in MeV
c     EGAM is gamma-ray energy in MeV
c
c   E1:  NOPTE1= 0: Single-particle approximation
c                1: Classical Lorentzian GDER
c                2: GDER with an energy and temperature dependent
c                   damping width (J.Kopecky, R.Chrien, Nucl.Phys.
c                   A468,p.285)
c                3: correct EGLO model - see 6
c                4: Kadmenskij-Markushev-Furman original Strength function
c                   (no high energy approximation according to Chrien)
c                5: The Chrien's Strength function (Nucl. Phys. A468, 285
c                   (1987)) only. In 3: is this model used only for
c                   high energy region
c                6: The strength function according to Chrien with
c                   phenomenological temperature dependent damping
c                   proposed by Kopecky (Distribution of Radiative Strength
c                   in Gd-156, 157 and 158 Nuclei)
c                   - I have found an error, correct EGLO is 3:
c                7: GDER with phenomenological temperature dependent
c                   damping width proposed by Kopecky
c                8: 4: with the first resonance of Lorentz type
c                9: 6: with the first resonance of Lorentz type
c               10: KMF (4:) for EG<4 MeV; lin. combination of KMF and BA
c                   for 4 MeV<EG<8 MeV; BA (1:) for EG>8 MeV
c               12: MLO1 (agrees with RIPL2 code)
c               16: MLO2 (agrees with RIPL2 code)
c               17: MLO3 (agrees with RIPL2 code)
c               31-40: correspond to 1-10 for high EGAM; for low EGAM original
c                      values are multiplied by a factor (given in input data)
c                      in between the PSF is a linear combination ... 
c                      motivation comes from Au
c               51: KMF (4:) without temperature-dependent term in damping
c                   width
c               52: EGLO (6:) without temperature-dependent term in damping
c                   width; temperature is taken into account only in the
c                   "second term" - FK*...
c               53: EELO (7:) without temperature-dependent term in damping
c                   width - i.e. no termperature dependence assumed
c               41: KMF according to Oslo group 
c
c   M1:  NOPTM1= 0: Single-particle approximation
c                1: Classical lorentzian GDMR
c                3: Scissors (first) resonance is build up only on states
c                   with excitation energy lower than PAR_M1(1)
c                4: Classical lorentzian build on the "background"
c                   that is described by the SP (constant function)
c               12: Scissors (first) resonance is build up only on states
c                   with excitation energy lower than PAR_M1(1) + SP on all states (originally NOPTM1=11)
c
c                ?: Enery of scissors resonance depends linearly on
c                   the energy of final state (and is build up only on states
c                   below certain excitation energy)
c                ?: Scissors resonance is considered only for primary transitions
c
c                ?: power dependence
c
c   E2:  NOPTE2= 0: Single-particle approximation
c                1: Classical Lorentzian GQER
c
C***********************************************************************
c
      PARAMETER  (PIH=  8.673592583E-08, ! 1/(3*(pi*hbar*c)**2)
     &            PIHQ= 5.204155555E-08, ! 1/(5*(pi*hbar*c)**2)
     &            PI42=39.4784176)       ! 4*pi**2
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000,
     *              MAXTPSF = 50     )   
      LOGICAL      LWARN_T
      REAL*4       GET_T_PSF2D
c
      COMMON /PHOTO/  NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     &                SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     &                DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     &                PAR_E1(3),DIPSLP,DIPZER
     &       /DEN/    NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *                DENLO,DENHI,DENPAR(4),
     *                ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
      COMMON /MAINE/  eall,EIN,efi,NOPTFL,BN,delta,nbin,ntotal,iregi,
     &                noptcs,nlinc,capfr(2)
C$OMP THREADPRIVATE(/MAINE/)
      COMMON /SFCE/   SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
     &                GE1DIS(2),GM1DIS(2),GE2DIS(2)
C$OMP THREADPRIVATE(/SFCE/)
      COMMON /incr/  nfilev,STEPS,enrgfin(18),buff(251),
     &               ncum,cumwidth(14),edo(18),eup(18)
      COMMON /PSF2DTAB/ TABTPSF(3,MAXTPSF),TABPSF2D(3,0:400,MAXTPSF),
     *                  TABPSF2D_LOG(3,0:400,MAXTPSF),
     *                  PSF2D_COEF(3,0:396,4,MAXTPSF),
     *                  NTPSF(3),LWARN_T(3)
c
      SGAMMA=0.
      SFCEE1=0
      SFCEM1=0
      SFCEE2=0
      IF ((ITYP.GT.4).OR.(ITYP.LT.1).OR.(EGAM.LE.0.)) RETURN
c
c
c*****                       E1 component
c
      IF (ITYP.EQ.1) THEN
c
        IF     (NOPTE1.EQ.0) THEN   ! The single-particle approximation
          SGAMMA=DEG*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.1) THEN   ! Classical Lorentzian
          Q=0.
          DO I=1,NGIGE              ! loop over both GDR peaks
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.2) THEN   ! GDER with E,T-dependent damping
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.3) THEN  !Empirical generalization of temperature
c          dependent damping according to Kopecky in Chrien model (EGLO)
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            WPHENZ=EK0-(1.-EK0)*EGZERO/(ER(I)-EGZERO)
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5*WPHENZ
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.4) THEN   ! Pure Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)          ! (no high energy approximation)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.5) THEN    !Pure Chrien model (similar to !OPT=3,
c          but no Kadmenskij for low EGAM) viz Kopecky
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.6) THEN  !Empirical generalization of temperature
c          dependent damping according to Kopecky in Chrien model (EGLO)
c          - an error found in the SLIM - corrected in NOPTE1=3
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.7) THEN   ! SMLO
          IF (EGAM.LE.EINI) THEN
            TFIN=SQRT((EINI-EGAM)/AMASS*10.0)
          ELSE
            TFIN=0.0
          ENDIF
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM*ER(I)+PI42*TFIN**2)/ER(I)**2
            QQ=SIG(I)*W0(I)*W*EGAM / (1.0-EXP(-EGAM/TFIN)) /
     &         ((EGAM**2-ER(I)**2)**2+(EGAM*W)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.207) THEN  !Empirical generelization of temperature
c        dependent damping according to Kopecky aplied to TD model (EELO)
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.8) THEN   ! Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)      ! with 1st resonance of Lorentz. shape
          Q=0.
          QQ=SIG(1)*
     &       (EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2))
          Q=Q+QQ
          DO I=2,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.9) THEN  !Empirical generalization of temperature
c          dependent damping according to Kopecky in Chrien model (EGLO)
c          with the first resonance of Lorentzian type
          TFIN=TERM(EINI-EGAM)
          Q=0.
          QQ=SIG(1)*
     &       (EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2))
          Q=Q+QQ
          DO I=2,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.10) THEN   ! KMF for low energies
          TFIN=TERM(EINI-EGAM)       ! Mix KMF and BA for higher energies
          Q=0.
          IF (EGAM.LE.PAR_E1(1)) THEN
           DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
           ENDDO
          ELSE
           x=(EGAM-PAR_E1(1))/(PAR_E1(2)-PAR_E1(1))         ! Admixture of BA to KMF
           IF (x.GT.1.) x=1.
           DO I=1,NGIGE
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+x*QQ
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+(1.-x)*QQ
           ENDDO
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.11) THEN                            ! Goriely 1D
          SGAMMA = APSF(EGAM,1)*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.12) THEN
C         2D PSF(Eg,Tf) or legacy 1D+enh
          IF (NTPSF(1).GE.2) THEN
            TFIN = GET_T_PSF2D(EINI-EGAM)
            SGAMMA = APSF2D(EGAM,TFIN,1)*EGAM**3
          ELSE
            SGAMMA = APSF(EGAM,1)
            SGAMMA = SGAMMA + EINI*PAR_E1(1)/(1.0+EXP(EGAM-PAR_E1(2)))
            SGAMMA = SGAMMA * EGAM**3
          ENDIF
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.25) THEN    !Pure Chrien+Kopecky model but different T (similar to OPT=3)
          TFIN=TERM(BN-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.31) THEN  ! Classical Lor.; suppressed for small EGAM
          Q=0.
          DO I=1,NGIGE              ! loop over both GDR peaks
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.34) THEN   ! Pure Fermi liquid theory (Kadmenskij)
c                                      suppressed for small EGAM
          TFIN=TERM(EINI-EGAM)          
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.36) THEN  !EGLO (6 - incorrect); suppressed for small EGAM
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.37) THEN  !EELO; suppressed for small EGAM
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.74) THEN   ! KMF with constant T
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.76) THEN  ! EGLO(6) with constant T
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.80) THEN   ! KMF for low energies
          TFIN=PAR_E1(3)                ! Mix KMF and BA for higher energies
          Q=0.
          IF (EGAM.LE.PAR_E1(1)) THEN
           DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
           ENDDO
          ELSE
           x=(EGAM-PAR_E1(1))/(PAR_E1(2)-PAR_E1(1))         ! Admixture of BA to KMF
           IF (x.GT.1.) x=1.
           DO I=1,NGIGE
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+x*QQ
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+(1.-x)*QQ
           ENDDO
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.83) THEN  ! EGLO(6) with constant T + Lorentzian pygmy (Oslo - actinides); PAR_E1(1) !!!
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=2,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          DO I=1,1              ! pygmy
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.82) THEN  ! EGLO(6) with constant T + Lorentzian pygmy (Oslo - actinides); PAR_E1(1) !!!
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=4,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          DO I=1,3              ! pygmy
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.84) THEN   ! GLO + gaussian PR (Oslo - Sn)
          Q=0.
          TFIN = PAR_E1(1)
          DO I = 2, NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          DO I = 1, 1
            QQ = SIG(I) / W0(I) / SQRT(2.0 * 3.141592) *
     &           EXP(-(EGAM-ER(I))**2 / 2.0 / W0(I)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA= PAR_E1(3) * PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

c
        ELSEIF (NOPTE1.EQ.64) THEN   ! SMLO (RIPL-3)
c
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W = EINI * W0(I) / ER(I)
            IF (TFIN.GT.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.67) THEN   ! Mughabghab+Dunford - formulation from RIPL2 
c
          S2=217.16/AMASS/AMASS      !global parametrization, see RIPL2 
          CDQ=1.05
          TFIN=TERM(EINI-EGAM)
          BETA2 = PAR_E1(3)*PAR_E1(3)    !square of deformation
c
          Q=0.
          DO I=1,NGIGE
            WDQ0= CDQ*SQRT(BETA2*ER(I)*ER(I)+ER(I)*S2)
            WDQ = CDQ*SQRT(BETA2*EGAM*EGAM+EGAM*S2)
            WC  = (W0(I)-WDQ0)/ER(I)/ER(I)*(EGAM**2+PI42*TFIN**2)
            W   = WDQ + WC 
            QQ=SIG(I)*W0(I)*FERMC*ER(I)*W/
     &         ((EGAM**2-ER(I)**2)**2+FERMC*(EGAM*W)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.68) THEN   ! Mughabghab+Dunford - original formulation (PLB487,155)  
c
          S2=217.16/AMASS/AMASS      !global parametrization, see RIPL2 
          CDQ=1.05
          TFIN=TERM(EINI-EGAM)
          BETA2 = PAR_E1(3)*PAR_E1(3)    !square of deformation
c
          Q=0.
          DO I=1,NGIGE
            WDQ0= CDQ*SQRT(BETA2*ER(I)*ER(I)+ER(I)*S2)
            WDQ = CDQ*SQRT(BETA2*EGAM*EGAM+EGAM*S2)
            WM  = (W0(I)-WDQ0)/ER(I)/ER(I)*(EGAM**2+PI42*TFIN**2)
            W   = WDQ + WM 
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.69) THEN   ! Goriely (PLB436,10) - in RIPL2
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*FERMC*(EGAM**2+PI42*TFIN**2)/ER(I)/EGAM
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+EGAM**2*W0(I)*W))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
c
        ELSEIF (NOPTE1.EQ.112) THEN   ! MLO1 (Original Plujko)
c         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
c         !ktere mohou podle 'Plujkovy teorie' nastat 
c         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
c         !hodnoty nekterych 'promennych' nastaveny zde
c
          TFIN=TERM(EINI-EGAM)
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          ER0PL=SMFREQ()
c
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)/ALPPL
            WDR=ER(I)**2/ALPPL
            WR=2*WDR*(ER(I)**2+ER0PL**2)/
     &         ((ER(I)**2-ER0PL**2)**2+4*(WDR*ER(I))**2)
            W=2*WD*W0(i)/WR*(ER(I)**2+ER0PL**2)/    
     &        ((ER(I)**2-ER0PL**2)**2+4*(WD*EGAM)**2)

            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.13) THEN   ! MLO1 (Original Plujko); Jina T!!!!
c         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
c         !ktere mohou podle 'Plujkovy teorie' nastat 
c         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
c         !hodnoty nekterych 'promennych' nastaveny zde
c
          TFIN=TERMDILG(EINI-EGAM)
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          ER0PL=SMFREQ()
c
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)/ALPPL
            WDR=ER(I)**2/ALPPL
            WR=2*WDR*(ER(I)**2+ER0PL**2)/
     &         ((ER(I)**2-ER0PL**2)**2+4*(WDR*ER(I))**2)
            W=2*WD*W0(i)/WR*(ER(I)**2+ER0PL**2)/    
     &        ((ER(I)**2-ER0PL**2)**2+4*(WD*EGAM)**2)

            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.14) THEN   ! MLO - RIPL2 
c         !it is not clear if there is a factor EINI or (EGAM+a*TFIN**2) 
c          or (EGAM+a*TINI**2) in the width WD 
c
          TFIN=TERM(EINI-EGAM)
          TINI=TERM(EINI)
c          TINI=TERM1(EINI)
          ALPPL=0.00542              !in RIPL2 error in the order (5.42e-4)
          FACTOR=1.0                 !All parameters recommended by Plujko
          ALPPL=ALPPL/FACTOR
          FKs0=0.3                          
          FNS=1.0                           
c
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)*ALPPL
c            WD=(EGAM+UINI)*ER(I)*ALPPL
            FKR=(W0(I)-ALPPL*ER(I)**2)/WWALL1()
            IF (EGAM.LE.(2*ER(I))) THEN
             FKS1=FKR+(FKs0-FKR)*ABS((EGAM-W0(I))/W0(I))**FNS             
            ELSE
             FKS1=FKs0
            ENDIF
            W=WD+FKS1*WWALL1()
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.16) THEN   ! MLO2 (Original Plujko)
c         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
c         !ktere mohou podle 'Plujkovy teorie' nastat 
c         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
c         !hodnoty nekterych 'promennych' nastaveny zde
c
          TFIN=TERM(EINI-EGAM)
          BC=1.0                            ! Recommended by Plujko
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          EFERMI=37.                        ! Recommended by Plujko 
          WWALL=6.857764*sqrt(EFERMI)/R0()  ! origin of the const. is a puzzle for me
          FKs0=0.3                          ! Recommended by Plujko
          FNS=1.0                           ! Recommended by Plujko
c
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)/ALPPL
            IF (NGIGE.EQ.1) THEN
              W=WD+FKS(EGAM,FNS,FKs0,WWALL,ALPPL)*WWALL
            ELSE
              W=WD+(W0(I)-ER(I)**2/ALPPL)*FKS(EGAM,FNS,FKs0,WWALL,ALPPL)
     &          /FKS(ER(I),FNS,FKs0,WWALL,ALPPL)
            ENDIF
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.17) THEN   ! MLO3 (Original Plujko)
c         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
c         !ktere mohou podle 'Plujkovy teorie' nastat 
c         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
c         !hodnoty nekterych 'promennych' nastaveny zde
c
          TFIN=TERM(EINI-EGAM)
          BC=1.0                            ! Recommended by Plujko
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          EFERMI=37.                        ! Recommended by Plujko 
          WWALL=6.857764*sqrt(EFERMI)/R0()  ! origin of the const. is a puzzle for me
          FKs0=0.3                          ! Recommended by Plujko
          FNS=1.0                           ! Recommended by Plujko
c
          Q=0.
          DO I=1,NGIGE
            WD=(EGAM**2+PI42*TFIN**2)/ALPPL
            IF (NGIGE.EQ.1) THEN
              W=WD+FKS(EGAM,FNS,FKs0,WWALL,ALPPL)*WWALL
            ELSE
              W=WD+(W0(I)-ER(I)**2/ALPPL)*FKS(EGAM,FNS,FKs0,WWALL,ALPPL)
     &          /FKS(ER(I),FNS,FKs0,WWALL,ALPPL)
            ENDIF
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.18) THEN   ! Sirotkin (1) - with KMF damping
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            O=FERMC*ER(I)
            IF (TFIN.NE.0) THEN
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/
     &         O**2/(1-exp(-EGAM/TFIN))
            ELSE
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/
     &         O**2
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.19) THEN   ! Sirotkin (2) - with damping from Sir
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
c            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            W=0.666*4./35.*(EGAM**2+PI42*TFIN**2)
            O=FERMC*ER(I)
            IF (TFIN.NE.0) THEN
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/
     &         O**2/(1-exp(-EGAM/TFIN))
            ELSE
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/
     &         O**2
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.20) THEN   ! Pure Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)       ! with theoretical damping from Sir
          Q=0.
          DO I=1,NGIGE
c            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            W=0.666*4./35.*(EGAM**2+PI42*TFIN**2)
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.51) THEN    ! Fermi liquid theory (Kadmenskij)
                                      ! without dependence on T !!!
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.52) THEN  !Modified empirical generalization
c          of temperature dependent damping according to Kopecky
c          in Chrien model (EGLO); Temperature not assumed in damping
c          width - see above
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2)/ER(I)**2
c          !     ^ .......... energy dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.53) THEN  !Modified empirical generelization
c        of temperature dependent damping according to Kopecky aplied
c        to TD model (EELO); temperature dependence removed - see above
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2)/ER(I)**2
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
c
c
        ELSEIF (NOPTE1.EQ.41) THEN   !!!!!! Oslo KMF for Yb
          TFIN=TERMOSLO(EINI-EGAM)        
c          TFIN=0.34        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.42) THEN   !!!!!! Oslo KMF - our temperature
          TFIN=TERM(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.43) THEN   !!!!!! Oslo KMF; including softpole
          TFIN=TERMOSLO(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.44) THEN   !!!!!! Oslo KMF - our temperature; including softpole
          TFIN=TERM(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.45) THEN   !!!!!! Oslo KMF - const T 
          TFIN=TERMOSLO(BN-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.46) THEN   !!!!!! Oslo KMF; including softpole
          TFIN=TERMOSLO(BN-EGAM)       !Temperature independent PSF  
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.47) THEN   !!!!!! Oslo KMF; including softpole
          TFIN=TERMOSLO(EINI-EGAM)     !!!! const below Etr   
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.48) THEN   !!!!!! Oslo KMF - T independent; including softpole
          TFIN=TERMOSLO(BN-EGAM)     !!!! const below Etr   
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN

        ELSEIF (NOPTE1.EQ.49) THEN   !!!!!! Oslo KMF - T const; including softpole
          TFIN=EK0        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN



        ELSEIF (NOPTE1.EQ.71) THEN   ! 34 (Kadmenskij)
c                                      with enhancement at high energies
          TFIN=TERM(EINI-EGAM)          
          Q=0.

          IF (EGAM.LT.7.35) THEN
           DO I=1,NGIGE
             W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
             QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
             Q=Q+QQ
           ENDDO
           x=DIPSLP*EGAM+DIPZER 
           if (x.LT.PAR_E1(3)) x=PAR_E1(3)
           if (x.GT.1.0)    x=1.0
           if ((EGAM.GT.6.0).AND.(EGAM.LT.7.35)) 
     &         x=x+x*2*(EGAM-6.0)
           SGAMMA=PIH*Q*EGAM**3*x
          ELSE
           DO I=1,NGIGE              ! BA
             QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
             Q=Q+QQ
           ENDDO
           SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.78) THEN   !!!!!! Voinov Mo - no LLR
c          TFIN=0.015        
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.79) THEN   !!!!!! Voinov Mo LLR
c          TFIN=0.015        
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=1,1                   ! Lorentzian LLR
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ/PAR_E1(1)
          ENDDO
          DO I=2,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
c
        ELSEIF (NOPTE1.EQ.28) THEN  !Hokus - pokus
          SGAMMA=DEG*EGAM**8
          SFCEE1=SGAMMA
          RETURN
        ENDIF
      ENDIF
c
c
c*****                       M1 component
c
      IF ((ITYP.EQ.2).OR.( ITYP.EQ.3)) THEN
c
        IF (NOPTM1.EQ.0) THEN       ! The single-particle approximation
          SGAMMA=DMG*EGAM*EGAM*EGAM
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.1) THEN   ! Classical Lorentzian
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.3) THEN   ! Scissors mode is build up only on states
          Q=0.                      ! with Efinal lower than Etr
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/
     &        ((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.112) THEN   ! Scissors mode is build up only on states
          Q=0.                       ! with Efinal lower than Etr; SP everywhere ! Originally NOPTM1=11
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/
     &        ((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.4) THEN   ! Classical Lorentzian on backgroung
          Q=0.                      ! which is described by SP up to Etr
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          IF (EGAM.LE.PAR_M1(1)) THEN
            SGAMMA=(PIH*Q+DMG)*EGAM**3
          ELSE
            SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.7) THEN   ! SMLO + "upbend"
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q
          SGAMMA = SGAMMA + 
     *             PAR_M1(1) * EXP(-PAR_M1(2)*EGAM) *
     *             (1.0+PAR_M1(3)*EGAM*EGAM*EGAM)
          SGAMMA = SGAMMA * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN


        ELSEIF (NOPTM1.EQ.5) THEN  ! SP on states
          Q=0.                     ! with Efinal lower than Etr
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &         ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=EINI-EGAM
          IF (Efinal.LE.PAR_M1(1)) THEN
            SGAMMA=(PIH*Q+DMG)*EGAM**3
          ELSE
            SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.6) THEN   ! Classical Lorentzian on backgroung
          Q=0.                      ! which is described by exponencial fcion
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=(PIH*Q+DMG*exp(-EGAM/PAR_M1(1)))*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.18) THEN   ! Classical Lorentzian on a background
          Q=0.                       ! given by exponencial fcion of Eg and Ei
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=EINI-EGAM
          SGAMMA=(PIH*Q+DMG*exp(-EINI/PAR_E1(3))
     *           * exp(-EGAM/PAR_M1(1)))*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
c

        ELSEIF (NOPTM1.EQ.31) THEN  ! SP on states
          Q=0.                      ! with Eini lower than Etr
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &         ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          IF (EINI.LE.PAR_M1(1)) THEN
            SGAMMA=(PIH*Q+DMG)*EGAM**3
          ELSE
            SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.32) THEN   ! Classical Lorentzian ...
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=(PIH*Q+DMG*exp(-EINI/PAR_M1(1)))*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN


c
c
        ELSEIF (NOPTM1.EQ.8) THEN   ! SP + Double-bumped Scissors mode - higher bump only on states
          Q=0.                      ! with Efinal lower than Etr
          DO I=1,1
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          DO I=3,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(2)*EGAM*WM0(2)**2/
     &        ((EGAM**2-ERM(2)**2)**2+(EGAM*WM0(2))**2)
            Q=Q+QQ
          ENDIF
c          SGAMMA=PIH*Q*EGAM**3
          SGAMMA=(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
c
        ELSEIF (NOPTM1.EQ.9) THEN   ! Scissors mode is build up only on
          Q=0.                      ! some of the terminal states
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF ((Efinal.EQ.enrgfin(1)).OR.(Efinal.EQ.enrgfin(2)).OR.
c     &        (Efinal.EQ.enrgfin(3)).OR.(Efinal.EQ.enrgfin(4)).OR.
c     &        (Efinal.EQ.enrgfin(5)).OR.(Efinal.EQ.enrgfin(6)).OR.
     &        (Efinal.EQ.enrgfin(3))) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/
     &        ((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.10) THEN   ! Scissors mode up to Etr 
          Q=0.                       ! on SP backgroung (everywhere)
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/
     &        ((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
c
        ELSEIF (NOPTM1.EQ.11) THEN                            ! Goriely 1D
          SGAMMA = APSF(EGAM,2)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.12) THEN
C         2D PSF(Eg,Tf) or legacy 1D+upbend
          IF (NTPSF(2).GE.2) THEN
            TFIN = GET_T_PSF2D(EINI-EGAM)
            SGAMMA = APSF2D(EGAM,TFIN,2)*EGAM**3
          ELSE
            SGAMMA = APSF(EGAM,2)
            SGAMMA = SGAMMA + 
     *               PAR_M1(1) * EXP(-PAR_M1(2)*EGAM) *
     *               (1.0+PAR_M1(3)*EGAM*EGAM*EGAM)
            SGAMMA = SGAMMA * EGAM**3
          ENDIF
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.13) THEN                            ! Goriely (incorrect Eg dep)
          SGAMMA = APSF(EGAM,2)
          SGAMMA = SGAMMA + 
     *          PAR_M1(1)*EXP(-PAR_M1(2)*EGAM)*(1.0+PAR_M1(3)*EGAM*EGAM)
          SGAMMA = SGAMMA * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
c

        ELSEIF (NOPTM1.EQ.14) THEN   ! Width of Scissors resonance depends on Efin
          Efinal=Eini-Egam           ! => total strength increases
          WM0SC=WM0(1)+PAR_M1(1)*Efinal
          Q=0.
          DO I=1,NGIGM
            IF (I.EQ.1) THEN
              WM0_I=WM0SC
            ELSE
              WM0_I=WM0(I)
            ENDIF
            QQ=SIGM(I)*EGAM*WM0_I**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0_I)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.15) THEN   ! Width of Scissors resonance depends on Efin
          Efinal=Eini-Egam           ! Total strength should remain 
          WM0SC=WM0(1)+PAR_M1(1)*Efinal
          SIGMSC=WM0(1)*SIGM(1)/WM0SC
          Q=0.
          DO I=1,NGIGM
            IF (I.EQ.1) THEN
              WM0_I=WM0SC
              SIGM_I=SIGMSC
            ELSE
              WM0_I=WM0(I)
              SIGM_I=SIGM(I)
            ENDIF
            QQ=SIGM_I*EGAM*WM0_I**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0_I)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
c
        ELSEIF (NOPTM1.EQ.16) THEN   ! \sigma_0 of Scissors resonance depends on Efin
          Efinal=Eini-Egam            
          SIGMSC = SIGM(1) + PAR_M1(1) * Efinal
          IF (SIGMSC.LE.0.0) SIGMSC = 0.0
          Q=0.
          DO I=1,NGIGM
            IF (I.EQ.1) THEN
              SIGM_I=SIGMSC
            ELSE
              SIGM_I=SIGM(I)
            ENDIF
            QQ=SIGM_I*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.21) THEN   ! Oslo - factor k
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.22) THEN   ! Oslo - factor k + softpole
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Q=Q+PAR_M1(2)/(EGAM**PAR_M1(3))
          SGAMMA=PAR_M1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.23) THEN   ! Oslo - factor k + softpole with a const below Etr
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.24) THEN   ! SP - factor k + softpole 
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.25) THEN   ! Oslo - factor k + softpole/rossendorf
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          QQ= PAR_M1(1) * EXP(-EGAM / PAR_E1(2)) 
          SGAMMA=PAR_E1(1) * (PIH*Q + QQ) * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.26) THEN   ! Oslo - factor k + softpole/rossendorf (no PAR_E1(1))
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          QQ= PAR_M1(1) * EXP(-EGAM / PAR_E1(2)) 
          SGAMMA= (PIH*Q + QQ) * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN


        ELSEIF (NOPTM1.EQ.27) THEN   ! Oslo - factor k not for LLR
          Q=0.
          DO I=1,1
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ/PAR_E1(1)
          ENDDO
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ENDIF
      ENDIF
c
c
c*****                       E2 component
c
      IF ((ITYP.EQ.2).OR.(ITYP.EQ.4)) THEN
        IF (NOPTE2.EQ.0) THEN
          SGAMMA=SGAMMA+QEL*EGAM**5
            SFCEE2=SGAMMA-SFCEM1
        ELSEIF (NOPTE2.EQ.11) THEN                            ! Goriely 1D
          SGAMMA = SGAMMA + APSF(EGAM,3)*EGAM**5
          SFCEE2=SGAMMA-SFCEM1
          RETURN
        ELSEIF (NOPTE2.EQ.12) THEN
C         2D PSF(Eg,Tf)
          IF (NTPSF(3).GE.2) THEN
            TFIN = GET_T_PSF2D(EINI-EGAM)
            SGAMMA = SGAMMA + APSF2D(EGAM,TFIN,3)*EGAM**5
          ELSE
            SGAMMA = SGAMMA + APSF(EGAM,3)*EGAM**5
          ENDIF
          SFCEE2=SGAMMA-SFCEM1
          RETURN
        ELSEIF ((NOPTE2.GE.13).AND.(NOPTE2.LE.15)) THEN
          WRITE(*,*) 'FATAL: Unsupported NOPTE2 = ', NOPTE2
          STOP 'FATAL: Unsupported NOPTE2'
        ELSEIF (NOPTE2.EQ.1) THEN                            ! Goriely
          Q=0.
          DO I=1,NGIGE2
            QQ=(SIGE(I)*WE0(I)**2)/(EGAM*((EGAM**2-ERE(I)**2)**2+
     &        (EGAM*WE0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=SGAMMA+PIHQ*Q*EGAM**5
          SFCEE2=SGAMMA-SFCEM1
        ENDIF
        RETURN
      ENDIF
c
      END
c
C***********************************************************************
      FUNCTION TERM(EEXC)
c
c     The effective temperature corresponding to the exc.energy EEXC
c
C***********************************************************************
C attention!! : PAIRING (not EONE) is used for the shift 
c
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
c
      EEFF=EEXC-PAIRING
      IF (EEFF.LT.0.) EEFF=0.
      TERM=SQRT((EEFF)/ASHELL)
      RETURN
      END
C
C***********************************************************************
      FUNCTION TERM1(EEXC)
c
c     The effective temperature corresponding to the exc.energy EEXC
c
C***********************************************************************
c
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
c
      EEFF=EEXC
      IF (EEFF.LT.0.) EEFF=0.
      TERM1=SQRT((EEFF)/ASHELL)
      RETURN
      END
C
C***********************************************************************
      FUNCTION TERMOSLO(EEXC)
c
c     The effective temperature corresponding to the exc.energy EEXC
c
C***********************************************************************
C attention!! : PAIRING (not EONE) used for shift
c
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
c
      EEFF=EEXC-PAIRING+6.6/(AMASS**0.32)
      IF (EEFF.LT.0.) EEFF=0.
      TERMOSLO=SQRT((EEFF)/ASHELL)
      RETURN
      END
C
C***********************************************************************
      FUNCTION EFEXC(EEXC)
c
c     The effective temperature corresponding to the exc.energy EEXC
c
C***********************************************************************
C attention!! : PAIRING (not EONE) used for shift
c
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
c
      EFEXC=EEXC-PAIRING
      RETURN
      END
C
C***********************************************************************
      FUNCTION TERMDILG(EEXC)
c
c     The effective temperature corresponding to the exc.energy EEXC
c
C***********************************************************************
C attention!! : PAIRING (not EONE) used for shift
c
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
c
      EEFF=EEXC-PAIRING
      IF (EEFF.LT.0.) EEFF=0.
      TERMDILG=(1.0 + SQRT(1.0+4.0*ASHELL*EEFF))/2.0/ASHELL
      RETURN
      END
C
C***********************************************************************
      FUNCTION SMFREQ()
c
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
c
      SMFREQ=41./AMASS**(0.33333)
      RETURN
      END
C
C***********************************************************************
      FUNCTION EREL1()
c
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
c
      AMASS13=AMASS**(0.33333)
      EREL1=31.2/AMASS13+20.6/sqrt(AMASS13)
      RETURN
      END
C
C***********************************************************************
      FUNCTION WREL1()
c
      WREL1=0.026*EREL1()**1.91
      RETURN
      END
C
C***********************************************************************
      FUNCTION R0()
c
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
c
      R0=1.27*AMASS**(0.33333)
      RETURN
      END
C
C***********************************************************************
      FUNCTION FKS(EX,FNS,FKs0,WWALL,ALPPL)
c
      EREL=EREL1()
      if (EX.GT.(2*EREL)) then
        FKS=FKs0     
      else
        FACT=WREL1()-EREL**2/ALPPL
        FKS=FACT/WWALL+(FKs0-FACT/WWALL)*ABS((EX-EREL)/EREL)**FNS
      endif
      RETURN
      END
C
C***********************************************************************
      FUNCTION WWALL1()
c
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
c
          WWALL1=36.43*AMASS**(-0.3333333) 
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE BKT_BUILD(T, N, NBK, IBK, E0, SC, LOK)
C***********************************************************************
C     Builds a uniform bucket index over a sorted table T(1..N) so that
C     #{i : T(i) < X} is obtained in O(1) operations with exact local
C     refinement in NLT_BK. Sets LOK=.FALSE. if T is not non-decreasing.
C
      INTEGER N, NBK, IBK(0:NBK), I, J
      REAL T(*)
      REAL*8 E0, SC, EB
      LOGICAL LOK
      LOK = .TRUE.
      DO I = 2, N
        IF (T(I).LT.T(I-1)) LOK = .FALSE.
      ENDDO
      E0 = DBLE(T(1))
      IF (N.LT.2 .OR. T(N).LE.T(1)) THEN
        SC = 0.D0
        DO I = 0, NBK
          IBK(I) = 0
        ENDDO
        RETURN
      ENDIF
      SC = DBLE(NBK) / (DBLE(T(N)) - E0)
      J = 0
      DO I = 0, NBK
        EB = E0 + DBLE(I) / SC
        DO WHILE (J.LT.N)
          IF (DBLE(T(J+1)).LT.EB) THEN
            J = J + 1
          ELSE
            EXIT
          ENDIF
        ENDDO
        IBK(I) = J
      ENDDO
      RETURN
      END
C
C***********************************************************************
      INTEGER FUNCTION NLT_BK(X, T, N, NBK, IBK, E0, SC)
C***********************************************************************
C     Returns the exact count #{i=1..N : T(i) < X} for a monotone table
C     T(1..N) using the uniform bucket guess from BKT_BUILD.
C
      INTEGER N, NBK, IBK(0:NBK), IB, J
      REAL T(*), X
      REAL*8 E0, SC, D
      D = (DBLE(X) - E0) * SC
      IF (D.LT.0.D0) THEN
        IB = 0
      ELSEIF (D.GE.DBLE(NBK)) THEN
        IB = NBK
      ELSE
        IB = INT(D)
      ENDIF
      J = IBK(IB)
      DO WHILE (J.LT.N)
        IF (T(J+1).LT.X) THEN
          J = J + 1
        ELSE
          EXIT
        ENDIF
      ENDDO
      DO WHILE (J.GT.0)
        IF (T(J).GE.X) THEN
          J = J - 1
        ELSE
          EXIT
        ENDIF
      ENDDO
      NLT_BK = J
      RETURN
      END
C
C***********************************************************************
      FUNCTION APSF(EGX,MTYP)
C***********************************************************************
C
C     This subroutine provides cubic interpolation of values
C     of PSFs - based on the function AICC.
C                                                 Version from 10-NOV-09
C
      DIMENSION XX(4),YY(4),A(4)
            INTEGER      NLMIN_TAB
      COMMON /PSFTAB/TABENPSF(3,0:400),TABPSF(3,0:400),NPSF(3)
      COMMON /LOGTAB/ TABENPSF_LOG(3,0:400),TABPSF_LOG(3,0:400),
     *                PSF_COEF(3,0:396,4),
     *                NLMIN_TAB(0:49,0:1)
      COMMON /PSFPC/  PSFC1(4,0:396,3),PSFC2(4,0:396,50,3)
      COMMON /BKTAB/  PSF_ENC(400,3),PSF_TTC(50,3),
     *                IBK_PSF(0:512,3),IBK_TT(0:64,3),IBK_LD(0:512),
     *                BK_E0(3),BK_SC(3),BK_T0(3),BK_TS(3),
     *                BK_L0,BK_LS,LBK_LD
      REAL*8          BK_E0,BK_SC,BK_T0,BK_TS,BK_L0,BK_LS
      INTEGER         IBK_PSF,IBK_TT,IBK_LD,NLT_BK
      LOGICAL         LBK_LD
C
      APSF = 0.0
c
c    Low level densities treated in a special way
c
      IF (EGX.LT.TABENPSF(MTYP,1)) GOTO 11            ! f(XL) = 0
      IF (EGX.GT.TABENPSF(MTYP,NPSF(MTYP))) GOTO 11   ! f(XL) = 0
      IF (EGX.EQ.EGX) THEN
        K = NLT_BK(EGX, PSF_ENC(1,MTYP), NPSF(MTYP), 512,
     *             IBK_PSF(0,MTYP), BK_E0(MTYP), BK_SC(MTYP))
        K = MIN(MAX(K-2, 0), NPSF(MTYP)-4)
        GOTO 1
      ENDIF

      IF (EGX.LE.TABENPSF(MTYP,2)) THEN
         K=0
         ELSE
         IF (EGX.LE.TABENPSF(MTYP,NPSF(MTYP)-1)) THEN
            ILO=3
            IHI=NPSF(MTYP)-1
            DO WHILE (ILO.LT.IHI)
               IMID=(ILO+IHI)/2
               IF (EGX.LE.TABENPSF(MTYP,IMID)) THEN
                  IHI=IMID
               ELSE
                  ILO=IMID+1
               ENDIF
            ENDDO
            K=ILO-3
            GO TO 1
            ELSE
            K=NPSF(MTYP)-4
         ENDIF
      ENDIF
C
    1 eglog=log(egx)
      DO J=1,4
         XX(J)=TABENPSF_LOG(MTYP,K+J)
      ENDDO !J
      APSF=0.
      DO I=1,4
         A(I)=PSFC1(I,K,MTYP)
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)*(XX(J)-EGLOG)
         ENDDO !J
         APSF=APSF+A(I)
      ENDDO !I
C
      APSF=exp(APSF)

   11 CONTINUE
      RETURN
      END
C
C***********************************************************************
C***********************************************************************
      REAL*4 FUNCTION GET_T_PSF2D(EEXC)
C***********************************************************************
C
C     Temperature helper isolated strictly for 2D PSF (APSF2D).
C     Preserves analytical models in TERM/TERM1/TERMOSLO/TERMDILG.
C
      REAL*4 EEXC, EEFF, ATEMP, TERM, TERMDILG
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
C
      IF (NOPTDE.EQ.14) THEN
        GET_T_PSF2D = ATEMP(EEXC)
      ELSEIF ((NOPTDE.EQ.0).OR.(NOPTDE.EQ.3)) THEN
        GET_T_PSF2D = TEMPER
      ELSEIF (NOPTDE.EQ.8) THEN
        GET_T_PSF2D = TEMPER09
      ELSEIF (NOPTDE.EQ.9) THEN
        EEFF = EEXC - PAIRING09
        IF (EEFF.LT.0.0) EEFF = 0.0
        IF (ASHELL09.GT.0.0) THEN
          GET_T_PSF2D = SQRT(EEFF / ASHELL09)
        ELSE
          GET_T_PSF2D = 0.0
        ENDIF
      ELSEIF (NOPTDE.EQ.7) THEN
        GET_T_PSF2D = TERMDILG(EEXC)
      ELSE
        EEFF = EEXC - PAIRING
        IF (EEFF.LT.0.0) EEFF = 0.0
        IF (ASHELL.GT.0.0) THEN
          GET_T_PSF2D = SQRT(EEFF / ASHELL)
        ELSE
          GET_T_PSF2D = 0.0
        ENDIF
      ENDIF
      RETURN
      END
C
C***********************************************************************
      FUNCTION APSF2D(EGX,TFIN,MTYP)
C***********************************************************************
C
C     2D interpolation of PSF f(E_gamma, T_f):
C     - Cubic Lagrange interpolation in ln(E_gamma) along T columns
C     - Linear interpolation in T_f between T_IT and T_{IT+1}
C     - Boundary clamping to [T_1, T_NTPSF] with a one-time warning
C     - Automatic fallback to non-negative bilinear interpolation on 2x2
C       cell if stencil contains non-positive (<=1.0E-30) points
C
      PARAMETER (MAXTPSF = 50)
      LOGICAL LWARN_T, L_FALLBACK, LWARN_READ
      DIMENSION XX(4),A(4)
      INTEGER NLMIN_TAB, IE_LO, IE_HI, IT_A, IT_B
      REAL WE, DE, F1, F2, F3, F4, F_LO, F_HI, WT, TFIN_LOC
      COMMON /PSFTAB/TABENPSF(3,0:400),TABPSF(3,0:400),NPSF(3)
      COMMON /LOGTAB/ TABENPSF_LOG(3,0:400),TABPSF_LOG(3,0:400),
     *                PSF_COEF(3,0:396,4),
     *                NLMIN_TAB(0:49,0:1)
      COMMON /PSF2DTAB/ TABTPSF(3,MAXTPSF),TABPSF2D(3,0:400,MAXTPSF),
     *                  TABPSF2D_LOG(3,0:400,MAXTPSF),
     *                  PSF2D_COEF(3,0:396,4,MAXTPSF),
     *                  NTPSF(3),LWARN_T(3)
      COMMON /PSFPC/  PSFC1(4,0:396,3),PSFC2(4,0:396,50,3)
      COMMON /BKTAB/  PSF_ENC(400,3),PSF_TTC(50,3),
     *                IBK_PSF(0:512,3),IBK_TT(0:64,3),IBK_LD(0:512),
     *                BK_E0(3),BK_SC(3),BK_T0(3),BK_TS(3),
     *                BK_L0,BK_LS,LBK_LD
      REAL*8          BK_E0,BK_SC,BK_T0,BK_TS,BK_L0,BK_LS
      INTEGER         IBK_PSF,IBK_TT,IBK_LD,NLT_BK
      LOGICAL         LBK_LD
C
      IF (NTPSF(MTYP).LE.1) THEN
        APSF2D = APSF(EGX, MTYP)
        RETURN
      ENDIF
      APSF2D = 0.0
      IF (EGX.LT.TABENPSF(MTYP,1)) RETURN
      IF (EGX.GT.TABENPSF(MTYP,NPSF(MTYP))) RETURN
C
C     Defensive check against NaN in TFIN
      IF (TFIN .NE. TFIN) THEN
        TFIN_LOC = TABTPSF(MTYP, 1)
      ELSE
        TFIN_LOC = TFIN
      ENDIF
C
      IF (EGX.EQ.EGX) THEN
        NPS = NPSF(MTYP)
        NLT = NLT_BK(EGX, PSF_ENC(1,MTYP), NPS, 512,
     *               IBK_PSF(0,MTYP), BK_E0(MTYP), BK_SC(MTYP))
        K = MIN(MAX(NLT-2, 0), NPS-4)
        IE_LO = NLT
        IF (PSF_ENC(MIN(NLT+1,NPS),MTYP).EQ.EGX) IE_LO = NLT + 1
        IE_LO = MIN(IE_LO, NPS-1)
        GOTO 6003
      ENDIF
C
C     Determine 4-point energy stencil index K
      IF (EGX.LE.TABENPSF(MTYP,2)) THEN
        K = 0
      ELSEIF (EGX.LE.TABENPSF(MTYP,NPSF(MTYP)-1)) THEN
        ILO = 3
        IHI = NPSF(MTYP) - 1
        DO WHILE (ILO.LT.IHI)
          IMID = (ILO + IHI) / 2
          IF (EGX.LE.TABENPSF(MTYP,IMID)) THEN
            IHI = IMID
          ELSE
            ILO = IMID + 1
          ENDIF
        ENDDO
        K = ILO - 3
      ELSE
        K = NPSF(MTYP) - 4
      ENDIF
C
C     Determine 2-point bounding cell in energy [IE_LO, IE_LO+1]
      IE_LO = 1
      IE_HI = NPSF(MTYP)
      DO WHILE (IE_LO.LT.IE_HI-1)
        IMID = (IE_LO + IE_HI) / 2
        IF (EGX.GE.TABENPSF(MTYP,IMID)) THEN
          IE_LO = IMID
        ELSE
          IE_HI = IMID
        ENDIF
      ENDDO
 6003 CONTINUE
      DE = TABENPSF(MTYP,IE_LO+1) - TABENPSF(MTYP,IE_LO)
      IF (DE.GT.0.0) THEN
        WE = (EGX - TABENPSF(MTYP,IE_LO)) / DE
      ELSE
        WE = 0.0
      ENDIF
      WE = MAX(0.0, MIN(1.0, WE))
C
      EGLOG = log(EGX)
      DO J = 1, 4
        XX(J) = TABENPSF_LOG(MTYP, K+J)
      ENDDO
C
      NT_M = NTPSF(MTYP)
      IF (TFIN_LOC.LE.TABTPSF(MTYP,1)) THEN
        IF (TFIN_LOC.LT.TABTPSF(MTYP,1)-1.0E-6) THEN
C$OMP ATOMIC READ
          LWARN_READ = LWARN_T(MTYP)
          IF (.NOT. LWARN_READ) THEN
C$OMP CRITICAL(WARN_CRIT)
            IF (.NOT. LWARN_T(MTYP)) THEN
              LWARN_T(MTYP) = .TRUE.
              WRITE(*,*) 'WARNING: T_f=', TFIN_LOC,
     *                   ' below 2D PSF T_min=', TABTPSF(MTYP,1),
     *                   ' for MTYP=', MTYP, ' (clamped).'
            ENDIF
C$OMP END CRITICAL(WARN_CRIT)
          ENDIF
        ENDIF
        IT_A = 1
        IT_B = 1
        WT = 0.0
      ELSEIF (TFIN_LOC.GE.TABTPSF(MTYP,NT_M)) THEN
        IF (TFIN_LOC.GT.TABTPSF(MTYP,NT_M)+1.0E-6) THEN
C$OMP ATOMIC READ
          LWARN_READ = LWARN_T(MTYP)
          IF (.NOT. LWARN_READ) THEN
C$OMP CRITICAL(WARN_CRIT)
            IF (.NOT. LWARN_T(MTYP)) THEN
              LWARN_T(MTYP) = .TRUE.
              WRITE(*,*) 'WARNING: T_f=', TFIN_LOC,
     *                   ' above 2D PSF T_max=', TABTPSF(MTYP,NT_M),
     *                   ' for MTYP=', MTYP, ' (clamped).'
            ENDIF
C$OMP END CRITICAL(WARN_CRIT)
          ENDIF
        ENDIF
        IT_A = NT_M
        IT_B = NT_M
        WT = 0.0
      ELSE
        ILO = NLT_BK(TFIN_LOC, PSF_TTC(1,MTYP), NT_M, 64,
     *               IBK_TT(0,MTYP), BK_T0(MTYP), BK_TS(MTYP))
        IT_A = ILO
        IT_B = ILO + 1
        WT = (TFIN_LOC - TABTPSF(MTYP,IT_A)) /
     *       (TABTPSF(MTYP,IT_B) - TABTPSF(MTYP,IT_A))
        WT = MAX(0.0, MIN(1.0, WT))
      ENDIF
C
C     Check for non-positive points in 4-point energy stencil
      L_FALLBACK = .FALSE.
      DO I = 1, 4
        IF (TABPSF2D(MTYP, K+I, IT_A) .LE. 1.0E-30) L_FALLBACK = .TRUE.
        IF (TABPSF2D(MTYP, K+I, IT_B) .LE. 1.0E-30) L_FALLBACK = .TRUE.
      ENDDO
C
      IF (L_FALLBACK) THEN
        F1 = MAX(0.0, TABPSF2D(MTYP, IE_LO,   IT_A))
        F2 = MAX(0.0, TABPSF2D(MTYP, IE_LO+1, IT_A))
        F_LO = F1 + WE * (F2 - F1)
        IF (IT_A.EQ.IT_B) THEN
          APSF2D = MAX(0.0, F_LO)
          RETURN
        ENDIF
        F3 = MAX(0.0, TABPSF2D(MTYP, IE_LO,   IT_B))
        F4 = MAX(0.0, TABPSF2D(MTYP, IE_LO+1, IT_B))
        F_HI = F3 + WE * (F4 - F3)
        APSF2D = MAX(0.0, F_LO + WT * (F_HI - F_LO))
        RETURN
      ENDIF
C
C     Cubic Lagrange interpolation in ln(EG) along bracketing columns
      S_LO = 0.0
      DO I = 1, 4
        A(I) = PSFC2(I, K, IT_A, MTYP)
        DO J = 1, 4
          IF (I.NE.J) A(I) = A(I) * (XX(J) - EGLOG)
        ENDDO
        S_LO = S_LO + A(I)
      ENDDO
      F_LO = exp(S_LO)
      IF (IT_A.EQ.IT_B) THEN
        APSF2D = F_LO
        RETURN
      ENDIF
C
      S_HI = 0.0
      DO I = 1, 4
        A(I) = PSFC2(I, K, IT_B, MTYP)
        DO J = 1, 4
          IF (I.NE.J) A(I) = A(I) * (XX(J) - EGLOG)
        ENDDO
        S_HI = S_HI + A(I)
      ENDDO
      F_HI = exp(S_HI)
C
      APSF2D = F_LO + WT * (F_HI - F_LO)
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE READ_PSF_FILE(MTYP, NAME)
C***********************************************************************
C
C     Reads 1D (NOPT=11) or 2D temperature-dependent (NOPT=12) PSF file
C     and precomputes logarithms and 4-point Lagrange weights.
C
      PARAMETER (MAXTPSF = 50)
      CHARACTER*12 NAME
      CHARACTER*1024 LINE_BUF
      INTEGER MTYP, NPSF_L, NTPSF_L, IOS_L, IEN, IT, NLMIN_TAB
      INTEGER NOPT_CUR
      LOGICAL LWARN_T
      COMMON /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     *               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     *               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     *               PAR_E1(3),DIPSLP,DIPZER
      COMMON /PSFTAB/TABENPSF(3,0:400),TABPSF(3,0:400),NPSF(3)
      COMMON /LOGTAB/ TABENPSF_LOG(3,0:400),TABPSF_LOG(3,0:400),
     *                PSF_COEF(3,0:396,4),
     *                NLMIN_TAB(0:49,0:1)
      COMMON /PSF2DTAB/ TABTPSF(3,MAXTPSF),TABPSF2D(3,0:400,MAXTPSF),
     *                  TABPSF2D_LOG(3,0:400,MAXTPSF),
     *                  PSF2D_COEF(3,0:396,4,MAXTPSF),
     *                  NTPSF(3),LWARN_T(3)
C
      NOPT_CUR = 0
      IF (MTYP.EQ.1) NOPT_CUR = NOPTE1
      IF (MTYP.EQ.2) NOPT_CUR = NOPTM1
      IF (MTYP.EQ.3) NOPT_CUR = NOPTE2
C
      OPEN (UNIT=5, FILE=NAME, STATUS='OLD')
      READ(5, '(A)') LINE_BUF
    5 READ(5, '(A)', IOSTAT=IOS_L) LINE_BUF
      IF (IOS_L.NE.0) THEN
        WRITE(*,*) 'FATAL: Unexpected EOF reading header in ', NAME
        STOP 'FATAL: Invalid PSF file header'
      ENDIF
      LINE_BUF = ADJUSTL(LINE_BUF)
      IF (LINE_BUF(1:1).EQ.'*' .OR. LINE_BUF(1:1).EQ.'#' .OR.
     *    LINE_BUF(1:1).EQ.'!' .OR. LINE_BUF(1:1).EQ.'c' .OR.
     *    LINE_BUF(1:1).EQ.'C' .OR. LEN_TRIM(LINE_BUF).EQ.0) GOTO 5
      NTPSF_L = 1
      READ(LINE_BUF, *, IOSTAT=IOS_L) NPSF_L, NTPSF_L
      IF (IOS_L.NE.0) THEN
        READ(LINE_BUF, *, IOSTAT=IOS_L) NPSF_L
        IF (IOS_L.NE.0) THEN
          WRITE(*,*) 'FATAL: Cannot parse grid counts from line: ',
     *               TRIM(LINE_BUF), ' in file: ', NAME
          STOP 'FATAL: Invalid PSF file header'
        ENDIF
        NTPSF_L = 1
      ENDIF
      IF (NOPT_CUR.EQ.11 .OR. NOPT_CUR.EQ.13) THEN
        NTPSF_L = 1
      ENDIF
      IF (NPSF_L.LT.4 .OR. NPSF_L.GT.400 .OR. NTPSF_L.LT.1 .OR.
     *    NTPSF_L.GT.MAXTPSF) THEN
        WRITE(*,*) 'FATAL: Invalid PSF grid dimensions in ', NAME,
     *             ': NPSF=', NPSF_L, ' NTPSF=', NTPSF_L
        STOP 'FATAL: Invalid PSF grid dimensions'
      ENDIF
      NPSF(MTYP) = NPSF_L
      NTPSF(MTYP) = NTPSF_L
      LWARN_T(MTYP) = .FALSE.
C
      IF (NTPSF(MTYP).LE.1) THEN
        NTPSF(MTYP) = 1
        READ(5, '(A)', IOSTAT=IOS_L) LINE_BUF
        IF (IOS_L.NE.0) THEN
          WRITE(*,*) 'FATAL: Unexpected EOF reading line 3 in ', NAME
          STOP 'FATAL: Unexpected EOF in 1D PSF'
        ENDIF
    6   READ(5, '(A)', IOSTAT=IOS_L) LINE_BUF
        IF (IOS_L.NE.0) THEN
          WRITE(*,*) 'FATAL: Unexpected EOF in 1D PSF data: ', NAME
          STOP 'FATAL: Unexpected EOF in 1D PSF'
        ENDIF
        LINE_BUF = ADJUSTL(LINE_BUF)
        IF (LINE_BUF(1:1).EQ.'*' .OR. LINE_BUF(1:1).EQ.'#' .OR.
     *      LINE_BUF(1:1).EQ.'!' .OR. LINE_BUF(1:1).EQ.'c' .OR.
     *      LINE_BUF(1:1).EQ.'C' .OR. LEN_TRIM(LINE_BUF).EQ.0) GOTO 6
        READ(LINE_BUF, *) TABENPSF(MTYP,1), TABPSF(MTYP,1)
        IF (TABENPSF(MTYP,1).LE.0.0)
     *    STOP 'FATAL: Non-positive PSF energy'
        TABENPSF_LOG(MTYP,1) = log(TABENPSF(MTYP,1))
        TABPSF_LOG(MTYP,1) = log(MAX(TABPSF(MTYP,1), 1.0E-35))
        DO IEN = 2, NPSF(MTYP)
          READ(5, *) TABENPSF(MTYP,IEN), TABPSF(MTYP,IEN)
          IF (TABENPSF(MTYP,IEN).LE.TABENPSF(MTYP,IEN-1)) THEN
            WRITE(*,*) 'FATAL: Non-increasing PSF energy in ', NAME,
     *                 ': IEN=', IEN, ' E=', TABENPSF(MTYP,IEN)
            STOP 'FATAL: Non-increasing PSF energy'
          ENDIF
          TABENPSF_LOG(MTYP,IEN) = log(TABENPSF(MTYP,IEN))
          TABPSF_LOG(MTYP,IEN) = log(MAX(TABPSF(MTYP,IEN), 1.0E-35))
        ENDDO
        CLOSE(5)
      ELSE
        IF (NTPSF(MTYP).GT.MAXTPSF) THEN
          WRITE(*,*) 'ERROR: NTPSF > MAXTPSF in ', NAME
          STOP 'ERROR'
        ENDIF
   10   READ(5, '(A)', IOSTAT=IOS_L) LINE_BUF
        IF (IOS_L.NE.0) THEN
          WRITE(*,*) 'FATAL: Unexpected EOF reading T grid in ', NAME
          STOP 'FATAL: Unexpected EOF in 2D PSF'
        ENDIF
        LINE_BUF = ADJUSTL(LINE_BUF)
        IF (LINE_BUF(1:1).EQ.'*' .OR. LINE_BUF(1:1).EQ.'#' .OR.
     *      LINE_BUF(1:1).EQ.'!' .OR. LINE_BUF(1:1).EQ.'c' .OR.
     *      LINE_BUF(1:1).EQ.'C' .OR. LEN_TRIM(LINE_BUF).EQ.0) GOTO 10
        READ(LINE_BUF, *) (TABTPSF(MTYP,IT), IT=1,NTPSF(MTYP))
        DO IT = 2, NTPSF(MTYP)
          IF (TABTPSF(MTYP,IT).LE.TABTPSF(MTYP,IT-1)) THEN
            WRITE(*,*) 'FATAL: Non-increasing PSF T grid in ', NAME,
     *                 ': IT=', IT, ' T=', TABTPSF(MTYP,IT)
            STOP 'FATAL: Non-increasing PSF T grid'
          ENDIF
        ENDDO
   20   READ(5, '(A)', IOSTAT=IOS_L) LINE_BUF
        IF (IOS_L.NE.0) THEN
          WRITE(*,*) 'FATAL: Unexpected EOF reading 2D data in ', NAME
          STOP 'FATAL: Unexpected EOF in 2D PSF'
        ENDIF
        LINE_BUF = ADJUSTL(LINE_BUF)
        IF (LINE_BUF(1:1).EQ.'*' .OR. LINE_BUF(1:1).EQ.'#' .OR.
     *      LINE_BUF(1:1).EQ.'!' .OR. LINE_BUF(1:1).EQ.'c' .OR.
     *      LINE_BUF(1:1).EQ.'C' .OR. LEN_TRIM(LINE_BUF).EQ.0) GOTO 20
        READ(LINE_BUF, *) TABENPSF(MTYP,1),
     *                    (TABPSF2D(MTYP,1,IT), IT=1,NTPSF(MTYP))
        IF (TABENPSF(MTYP,1).LE.0.0)
     *    STOP 'FATAL: Non-positive PSF energy'
        DO IEN = 2, NPSF(MTYP)
          READ(5, *) TABENPSF(MTYP,IEN),
     *               (TABPSF2D(MTYP,IEN,IT), IT=1,NTPSF(MTYP))
          IF (TABENPSF(MTYP,IEN).LE.TABENPSF(MTYP,IEN-1)) THEN
            WRITE(*,*) 'FATAL: Non-increasing PSF energy in ', NAME,
     *                 ': IEN=', IEN, ' E=', TABENPSF(MTYP,IEN)
            STOP 'FATAL: Non-increasing PSF energy'
          ENDIF
        ENDDO
        CLOSE(5)
        DO IEN = 1, NPSF(MTYP)
          TABENPSF_LOG(MTYP,IEN) = log(TABENPSF(MTYP,IEN))
          TABPSF(MTYP,IEN) = TABPSF2D(MTYP,IEN,1)
          TABPSF_LOG(MTYP,IEN) = log(MAX(TABPSF(MTYP,IEN), 1.0E-35))
          DO IT = 1, NTPSF(MTYP)
            IF (TABPSF2D(MTYP,IEN,IT).GT.0.0) THEN
              TABPSF2D_LOG(MTYP,IEN,IT) = log(TABPSF2D(MTYP,IEN,IT))
            ELSE
              TABPSF2D_LOG(MTYP,IEN,IT) = log(1.0E-35)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      CALL PSF_PRECOMP(MTYP)
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE PSF_PRECOMP(MTYP)
C***********************************************************************
C     Precomputes uniform bucket indices and 4-point cubic Lagrange
C     weights in ln(E) for APSF (PSFC1) and APSF2D (PSFC2) using the
C     exact division order of the on-the-fly routines. Guards non-positive
C     1D entries (e.g. zero entries in column T_1=0 of a 2D PSF table)
C     against log(0) exceptions.
C
      PARAMETER (MAXTPSF = 50)
      INTEGER NLMIN_TAB
      LOGICAL LWARN_T, LOK
      DIMENSION XX(4),YY(4)
      COMMON /PSFTAB/TABENPSF(3,0:400),TABPSF(3,0:400),NPSF(3)
      COMMON /LOGTAB/ TABENPSF_LOG(3,0:400),TABPSF_LOG(3,0:400),
     *                PSF_COEF(3,0:396,4),
     *                NLMIN_TAB(0:49,0:1)
      COMMON /PSF2DTAB/ TABTPSF(3,MAXTPSF),TABPSF2D(3,0:400,MAXTPSF),
     *                  TABPSF2D_LOG(3,0:400,MAXTPSF),
     *                  PSF2D_COEF(3,0:396,4,MAXTPSF),
     *                  NTPSF(3),LWARN_T(3)
      COMMON /PSFPC/  PSFC1(4,0:396,3),PSFC2(4,0:396,50,3)
      COMMON /BKTAB/  PSF_ENC(400,3),PSF_TTC(50,3),
     *                IBK_PSF(0:512,3),IBK_TT(0:64,3),IBK_LD(0:512),
     *                BK_E0(3),BK_SC(3),BK_T0(3),BK_TS(3),
     *                BK_L0,BK_LS,LBK_LD
      REAL*8          BK_E0,BK_SC,BK_T0,BK_TS,BK_L0,BK_LS
      INTEGER         IBK_PSF,IBK_TT,IBK_LD
      LOGICAL         LBK_LD
C
      DO I = 1, NPSF(MTYP)
        PSF_ENC(I,MTYP) = TABENPSF(MTYP,I)
      ENDDO
      CALL BKT_BUILD(PSF_ENC(1,MTYP), NPSF(MTYP), 512,
     *               IBK_PSF(0,MTYP), BK_E0(MTYP), BK_SC(MTYP), LOK)
      DO I = 1, MAX(NTPSF(MTYP),1)
        PSF_TTC(I,MTYP) = TABTPSF(MTYP,I)
      ENDDO
      CALL BKT_BUILD(PSF_TTC(1,MTYP), MAX(NTPSF(MTYP),1), 64,
     *               IBK_TT(0,MTYP), BK_T0(MTYP), BK_TS(MTYP), LOK)
C
      DO K = 0, NPSF(MTYP)-4
        DO J = 1, 4
          XX(J) = log(TABENPSF(MTYP,K+J))
          IF (TABPSF(MTYP,K+J).GT.0.0) THEN
            YY(J) = log(TABPSF(MTYP,K+J))
          ELSE
            YY(J) = log(1.0E-35)
          ENDIF
        ENDDO
        DO I = 1, 4
          A = YY(I)
          DO J = 1, 4
            IF (I.NE.J) A = A/(XX(J)-XX(I))
          ENDDO
          PSFC1(I,K,MTYP) = A
        ENDDO
        IF (NTPSF(MTYP).GT.1) THEN
          DO J = 1, 4
            XX(J) = TABENPSF_LOG(MTYP,K+J)
          ENDDO
          DO IT = 1, NTPSF(MTYP)
            DO I = 1, 4
              A = TABPSF2D_LOG(MTYP,K+I,IT)
              DO J = 1, 4
                IF (I.NE.J) A = A/(XX(J)-XX(I))
              ENDDO
              PSFC2(I,K,IT,MTYP) = A
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
C
c***********************************************************************
      SUBROUTINE RANDINIT()
c***********************************************************************
c     Matrix of 4xNNR random seeds is generated 
      INTEGER*4  IRINIT,IDUMMY
      REAL       dummy

      COMMON /rnd/ IRINIT(4,999)
     &       /readi/ nevents,NSREAL,NRL,ecrit,xrayk,xrayl,numlev,
     &               elowlev(99),elowsp(99),ilowip(99),
     &               STDISa(0:2,0:30,-2:2,0:1),RADW

      DO I = 1, 4
        IDUMMY = IRINIT(I,1)
        DO J = 2, NSREAL
          DO K = 1, 100  
            dummy = RAN0(IDUMMY)
          ENDDO
          IRINIT(I,J) = IDUMMY - J
          IF (IRINIT(I,J).LE.0) IRINIT(I,J) = IRINIT(I,J) + 2 * J
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE
!***********************************************************************


