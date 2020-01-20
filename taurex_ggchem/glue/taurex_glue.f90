module fort_ggchem
  implicit none



  contains
      SUBROUTINE INIT_LEAN()
        use PARAMETERS,ONLY: elements
        use DUST_DATA,ONLY: NELEM,eps0,mass,muH,elnam,amu
        use EXCHANGE,ONLY: H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl, & 
                          Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge, & 
                           As,Se,Br,Kr,Rb,Sr,Y,Zr,W 
      implicit none
      integer,parameter :: qp=selected_real_kind(33,4931)
      integer :: i,j,nr
      real*8 :: abund(74,4),eps(NELEM),epsH,mfrac(NELEM)
      real*8 :: m,val,addH2O
      character(len=2) :: el
      character(len=20) :: elname
      character(len=10) :: source(4)
      character(len=200) :: line
      logical :: found
      
      elements     = 'H He C N O Na Mg Si Fe Al Ca Ti S Cl K Li el'

      write(*,*) 
      write(*,*) "elemental abundances and masses ..."
      write(*,*) "==================================="
      elnam(1)  = 'H '
      elnam(2)  = 'He'
      elnam(3)  = 'Li'
      elnam(4)  = 'Be'
      elnam(5)  = 'B '
      elnam(6)  = 'C '
      elnam(7)  = 'N '
      elnam(8)  = 'O '
      elnam(9)  = 'F '
      elnam(10) = 'Ne'
      elnam(11) = 'Na'
      elnam(12) = 'Mg'
      elnam(13) = 'Al'
      elnam(14) = 'Si'
      elnam(15) = 'P '
      elnam(16) = 'S '
      elnam(17) = 'Cl'
      elnam(18) = 'Ar'
      elnam(19) = 'K '
      elnam(20) = 'Ca'
      elnam(21) = 'Sc'
      elnam(22) = 'Ti'
      elnam(23) = 'V '
      elnam(24) = 'Cr'
      elnam(25) = 'Mn'
      elnam(26) = 'Fe'
      elnam(27) = 'Co'
      elnam(28) = 'Ni'
      elnam(29) = 'Cu'
      elnam(30) = 'Zn'
      elnam(31) = 'Ga'
      elnam(32) = 'Ge'
      elnam(33) = 'As'
      elnam(34) = 'Se'
      elnam(35) = 'Br'
      elnam(36) = 'Kr'
      elnam(37) = 'Rb'
      elnam(38) = 'Sr'
      elnam(39) = 'Y '
      elnam(40) = 'Zr'
      elnam(41) = 'W '

      mass(H)  = 1.008  * amu  
      mass(He) = 4.0026 * amu
      mass(Li) = 6.94   * amu  
      mass(Be) = 9.0122 * amu
      mass(B)  = 10.81  * amu  
      mass(C)  = 12.011 * amu  
      mass(N)  = 14.007 * amu  
      mass(O)  = 15.999 * amu  
      mass(F)  = 18.998 * amu
      mass(Ne) = 20.180 * amu 
      mass(Na) = 22.990 * amu
      mass(Mg) = 24.305 * amu 
      mass(Al) = 26.982 * amu
      mass(Si) = 28.085 * amu  
      mass(P)  = 30.974 * amu
      mass(S)  = 32.06  * amu  
      mass(Cl) = 35.45  * amu  
      mass(Ar) = 39.948 * amu  
      mass(K)  = 39.098 * amu 
      mass(Ca) = 40.078 * amu  
      mass(Sc) = 44.956 * amu
      mass(Ti) = 47.867 * amu  
      mass(V)  = 50.942 * amu 
      mass(Cr) = 51.996 * amu 
      mass(Mn) = 54.938 * amu
      mass(Fe) = 55.845 * amu  
      mass(Co) = 58.933 * amu
      mass(Ni) = 58.693 * amu 
      mass(Cu) = 63.546 * amu  
      mass(Zn) = 65.38  * amu  
      mass(Ga) = 69.723 * amu  
      mass(Ge) = 72.63  * amu  
      mass(As) = 74.922 * amu
      mass(Se) = 78.96  * amu  
      mass(Br) = 79.904 * amu  
      mass(Kr) = 83.798 * amu  
      mass(Rb) = 85.468 * amu 
      mass(Sr) = 87.62  * amu  
      mass(Y ) = 88.906 * amu
      mass(Zr) = 91.224 * amu  
      mass(W ) = 183.84 * amu       

      eps(H)  = 12.00 ! D0
      eps(He) = 10.99 ! D0
      eps(Li) =  1.16 ! D0
      eps(C)  =  8.55 ! D0
      eps(N)  =  7.97 ! D0
      eps(O)  =  8.87 ! D0
      eps(Ne) =  8.08 ! D0
      eps(Na) =  6.33 ! D0
      eps(Mg) =  7.58 ! D0
      eps(Al) =  6.47 ! D0
      eps(Si) =  7.55 ! D0
      eps(S)  =  7.33 ! D0
      eps(Cl) =  5.50 ! D0
      eps(K)  =  5.12 ! D0
      eps(Ca) =  6.36 ! D0
      eps(Ti) =  5.02 ! D0
      eps(Cr) =  5.67 ! D0
      eps(Mn) =  5.39 ! D0
      eps(Fe) =  7.50 ! D0
      eps(Ni) =  6.25 ! D0

      eps(H)  = 12.00    
      eps(He) = 10.93  
      eps(Li) = 1.05     
      eps(Be) = 1.38   
      eps(B)  = 2.70     
      eps(C)  = 8.43     
      eps(N)  = 7.83     
      eps(O)  = 8.69     
      eps(F)  = 4.56  
      eps(Ne) = 7.93    
      eps(Na) = 6.24  
      eps(Mg) = 7.60    
      eps(Al) = 6.45  
      eps(Si) = 7.51     
      eps(P)  = 5.41  
      eps(S)  = 7.12     
      eps(Cl) = 5.50     
      eps(Ar) = 6.40     
      eps(K)  = 5.03    
      eps(Ca) = 6.34     
      eps(Sc) = 3.15  
      eps(Ti) = 4.95     
      eps(V)  = 3.93    
      eps(Cr) = 5.64    
      eps(Mn) = 5.43  
      eps(Fe) = 7.50     
      eps(Co) = 4.99  
      eps(Ni) = 6.22    
      eps(Cu) = 4.19     
      eps(Zn) = 4.56     
      eps(Ga) = 3.04     
      eps(Ge) = 3.65     
      eps(As) = -40.   
      eps(Se) = -40.     
      eps(Br) = -40.     
      eps(Kr) = 3.25     
      eps(Rb) = 2.52    
      eps(Sr) = 2.87     
      eps(Y ) = 2.21  
      eps(Zr) = 2.58  
      eps(W ) = 0.85   

      do i=1,NELEM
        eps(i) = 10.Q0 ** (eps(i)-12.Q0)
      enddo
      eps0 = eps
      
      end subroutine



      SUBROUTINE INIT_TAUREX_CHEMISTRY(ndispol, dispolfiles, do_charge)

        use PARAMETERS,ONLY: elements,initchem_info
        use CHEMISTRY,ONLY: NMOLdim,NMOLE,NELM,catm,cmol,el,dispol_file, &
           source,fit,natom,a,error,i_nasa, &
            m_kind,m_anz,elnum,elion,charge, &
            el,H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc, &
            Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,W
        use DUST_DATA,ONLY: mass,mel,amu
        use EXCHANGE,ONLY: nmol,mmol

        implicit none

        character*200, intent(in) :: dispolfiles(ndispol)
        character(len=2) :: cel(40),elnam
        integer, intent(in) :: ndispol
        logical, intent(in) :: do_charge
        character(len=1024) :: filename
        character(len=1024) :: line
        character(len=20) :: molname,upper,leer=''
        logical :: found,charged
        real*8 :: fiterr


        integer :: loop,i,ii,j,iel,e,smax,ret


        charge = do_charge
        cel(:) = '.'
        read(elements,*,end=100) cel
   100  NELM = 0
        charge = .false.
        do i=1,99
          if (cel(i)=='.' .or. cel(i) == ' ') exit
          elnam = trim(cel(i))
          if (elnam == '') exit
          found = .false.
          do e=1,NELM
            if (elnam.eq.catm(e)) then
              found=.true.
              exit
            endif  
          enddo
          NELM = NELM+1
          catm(NELM) = elnam
          if     (elnam=='el') then; el=NELM ; charge=.true.
          elseif (elnam=='H')  then;  H=NELM ; elnum(NELM)=1 
          elseif (elnam=='He') then; He=NELM ; elnum(NELM)=2 
          elseif (elnam=='Li') then; Li=NELM ; elnum(NELM)=3
          elseif (elnam=='Be') then; Be=NELM ; elnum(NELM)=4 
          elseif (elnam=='B')  then;  B=NELM ; elnum(NELM)=5 
          elseif (elnam=='C')  then;  C=NELM ; elnum(NELM)=6 
          elseif (elnam=='N')  then;  N=NELM ; elnum(NELM)=7 
          elseif (elnam=='O')  then;  O=NELM ; elnum(NELM)=8 
          elseif (elnam=='F')  then;  F=NELM ; elnum(NELM)=9 
          elseif (elnam=='Ne') then; Ne=NELM ; elnum(NELM)=10
          elseif (elnam=='Na') then; Na=NELM ; elnum(NELM)=11 
          elseif (elnam=='Mg') then; Mg=NELM ; elnum(NELM)=12
          elseif (elnam=='Al') then; Al=NELM ; elnum(NELM)=13
          elseif (elnam=='Si') then; Si=NELM ; elnum(NELM)=14
          elseif (elnam=='P')  then;  P=NELM ; elnum(NELM)=15 
          elseif (elnam=='S')  then;  S=NELM ; elnum(NELM)=16
          elseif (elnam=='Cl') then; Cl=NELM ; elnum(NELM)=17
          elseif (elnam=='Ar') then; Ar=NELM ; elnum(NELM)=18
          elseif (elnam=='K')  then;  K=NELM ; elnum(NELM)=19
          elseif (elnam=='Ca') then; Ca=NELM ; elnum(NELM)=20
          elseif (elnam=='Sc') then; Sc=NELM ; elnum(NELM)=21
          elseif (elnam=='Ti') then; Ti=NELM ; elnum(NELM)=22
          elseif (elnam=='V')  then;  V=NELM ; elnum(NELM)=23
          elseif (elnam=='Cr') then; Cr=NELM ; elnum(NELM)=24
          elseif (elnam=='Mn') then; Mn=NELM ; elnum(NELM)=25
          elseif (elnam=='Fe') then; Fe=NELM ; elnum(NELM)=26
          elseif (elnam=='Co') then; Co=NELM ; elnum(NELM)=27
          elseif (elnam=='Ni') then; Ni=NELM ; elnum(NELM)=28
          elseif (elnam=='Cu') then; Cu=NELM ; elnum(NELM)=29
          elseif (elnam=='Zn') then; Zn=NELM ; elnum(NELM)=30
          elseif (elnam=='Ga') then; Ga=NELM ; elnum(NELM)=31
          elseif (elnam=='Ge') then; Ge=NELM ; elnum(NELM)=32 
          elseif (elnam=='As') then; As=NELM ; elnum(NELM)=33 
          elseif (elnam=='Se') then; Se=NELM ; elnum(NELM)=34 
          elseif (elnam=='Br') then; Br=NELM ; elnum(NELM)=35 
          elseif (elnam=='Kr') then; Kr=NELM ; elnum(NELM)=36 
          elseif (elnam=='Rb') then; Rb=NELM ; elnum(NELM)=37 
          elseif (elnam=='Sr') then; Sr=NELM ; elnum(NELM)=38 
          elseif (elnam=='Y')  then;  Y=NELM ; elnum(NELM)=39 
          elseif (elnam=='Zr') then; Zr=NELM ; elnum(NELM)=40
          elseif (elnam=='W')  then;  W=NELM ; elnum(NELM)=41
          else
            print*,'element '//elnam
            stop "*** unknown element "
          endif
          if (initchem_info) then
            print*,'element '//elnam,elnum(NELM)
          endif
        enddo

          NMOLdim = 10000
          if (allocated(cmol)) then 
            deallocate(cmol)
          endif
          if (allocated(fit)) then 
            deallocate(fit)
          endif
          if (allocated(natom)) then 
            deallocate(natom)
          endif
          if (allocated(error)) then 
            deallocate(error)
          endif
          if (allocated(source)) then 
            deallocate(source)
          endif
          if (allocated(m_kind)) then 
            deallocate(m_kind)
          endif
          if (allocated(m_anz)) then 
            deallocate(m_anz)
          endif 

          if (allocated(a)) then 
            deallocate(a)
          endif

          allocate(cmol(NMOLdim),fit(NMOLdim),natom(NMOLdim))
          allocate(error(NMOLdim),a(NMOLdim,0:13))
          allocate(source(NMOLdim),m_kind(0:6,NMOLdim),m_anz(6,NMOLdim))

          a = 0.0
          i=1
          i_nasa = 0
          do loop=1,ndispol
            filename = trim(dispolfiles(loop))
            if (filename=='') exit
            filename = trim(filename)
            !if (initchem_info) write(*,*)
            write(*,*) 'reading kp-data from ' &
                     //trim(filename)//" ..."
            open(unit=12, file=filename, status='old')
            read(12,*) NMOLdim
            do ii=1,NMOLdim
              read(12,'(A300)') line
              read(line,*) molname,iel,cel(1:iel),m_anz(1:iel,i)
              molname=trim(molname)
              fiterr = 0.0
              j = index(line,"+/-")
              if (j>0) read(line(j+3:),*) fiterr
              error(i) = fiterr
              read(12,'(A300)') line
              read(line,*) fit(i)
              if (fit(i)==6) then
                read(line,*) fit(i),(a(i,j),j=0,7)
              elseif(fit(i)==7) then
                 i_nasa = 1
                 read(line,*) fit(i),(a(i,j),j=0,13)
              else   
                read(line,*) fit(i),(a(i,j),j=0,4)
              endif  
              m_kind(0,i) = iel
              natom(i) = 0
              found = .true.
              smax  = 0
              do j=1,m_kind(0,i)
                natom(i) = natom(i)+m_anz(j,i)
                if (index(elements,cel(j))<=0) found=.false. 
                smax = MAX(smax,ABS(m_anz(j,i)))
              enddo  
              if (.not.found) cycle    ! molecule has unselected element 
              if (smax>16) cycle       ! stoichiometric coefficient > 16
              if (m_kind(0,i)==1.and.natom(i)==1) cycle  ! pure atom
              j = index(molname,"_")
              if (j>1) then
                cmol(i) = upper(molname(j+1:)//leer(1:j))
              else
                cmol(i) = upper(molname)
              endif
              charged = .false.
              do j=1,m_kind(0,i)
                elnam = cel(j)
                found = .false.
                do e=1,NELM
                  if (elnam.eq.catm(e)) then
                    found=.true.
                    exit
                  endif  
                enddo
                if (.not.found) stop "*** should not occur"
                m_kind(j,i) = e
                if (e==el) charged=.true. 
              enddo  
              if (fit(i)==6.and.charged) cycle ! old charged BarklemCollet
              source(i) = loop
              call CHECK_DOUBLE(ndispol,dispolfiles,cmol(i),m_kind(:,i),m_anz(:,i),i,loop,ret)
              if (ret>0) then
                source(ret) = loop
                cmol(ret) = cmol(i)
                fit(ret)  = fit(i)
                a(ret,:)  = a(i,:)
                error(ret)= error(i)
                if (initchem_info) then
                  write(line,'(I4,A20,1x,99(I3,1x,A2,1x))') &
                    ret,trim(cmol(ret)),(m_anz(j,ret),cel(j),j=1,iel)
                  print*,trim(line)//"    OVERWRITE"
                endif  
              else  
                if (initchem_info) then
                  write(line,'(I4,A20,1x,99(I3,1x,A2,1x))') &
                     i,trim(cmol(i)),(m_anz(j,i),catm(m_kind(j,i)),j=1,iel)
                  if (loop==1) then
                    print*,trim(line)
                  else
                    print*,trim(line)//"    --> NEW" 
                  endif
                endif  
                if (iel==2.and. &
                ((m_kind(1,i)==el.and.m_anz(1,i)==-1.and.m_anz(2,i)==1).or. &
                 (m_kind(2,i)==el.and.m_anz(2,i)==-1.and.m_anz(1,i)==1)) &
                 ) then
                  e = m_kind(1,i)
                  if (e==el) e=m_kind(2,i)
                  elion(e) = i
                endif
                i = i+1
              endif  
            enddo
     200    close(12)
          enddo  
          NMOLE = i-1
          if (allocated(nmol)) then
              deallocate(nmol)
          endif

          if (allocated(mmol)) then
            deallocate(mmol)
          endif

          allocate(nmol(NMOLE),mmol(NMOLE))
    
          if(i_nasa==1) call NASA_POLYNOMIAL !Added by Yui Kawashima
    
          if (loop>1.and.initchem_info) then
            print* 
            do i=1,NMOLE
              print*,i,cmol(i),' ->  '//trim(dispolfiles(source(i)))
            enddo
          endif  
      
          !open(unit=1,file='chemicals.tex')
          !write(1,*) NMOLE
          !do i=1,NMOLE
          !  if (error(i)>0.0) then 
          !    write(1,3000)
         !      i,cmol(i),source(i),fit(i),a(i,0:4),error(i)
          !  else  
          !    write(1,3010)
         !      i,cmol(i),source(i),fit(i),a(i,0:4)
          !  endif  
          !enddo  
          !close(1)
          !stop
    
          do i=1,NMOLE
            mmol(i) = 0.d0
            do j=1,m_kind(0,i)
              if (m_kind(j,i)==el) then
                mmol(i) = mmol(i) + m_anz(j,i)*mel
              else
                mmol(i) = mmol(i) + m_anz(j,i)*mass(elnum(m_kind(j,i)))
              endif
            enddo
            !print*,cmol(i),mmol(i)/amu
          enddo  
    
          print* 
          print*,NMOLE,' species'
          print*,NELM,' elements'
          print'(99(A4))',(trim(catm(j)),j=1,NELM)
          print'(99(I4))',elnum(1:NELM)
          !print'(99(I4))',H,He,C,N,O,Si,Mg,Fe,Na,Al,S,Ca,Ti,Cl,K,Li,el
          if (charge) then
            print'(1x,99(A4))',(trim(cmol(elion(j))),j=1,el-1),'  ', &
                              (trim(cmol(elion(j))),j=el+1,NELM)
          endif  
          
     3000 format(I4," & ",A12," & (",I1,") & ",I1," & ", &
                5(1pE12.5," & "),"$\pm$",0pF4.2,"\\")
     3010 format(I4," & ",A12," & (",I1,") & ",I1," & ", &
                5(1pE12.5," & "),"\\")
          end
    
    !************************************************************************
          subroutine CHECK_DOUBLE(ndispol, dispolfiles,molname,kind,anz,N,loop,ret)
    !************************************************************************
          use PARAMETERS,ONLY: initchem_info
          use CHEMISTRY,ONLY: cmol,m_kind,m_anz,source
          implicit none
          character(len=*), intent(in) :: dispolfiles(ndispol)
          character(len=20) :: molname
          integer,intent(IN) :: kind(0:6),anz(6),N,loop, ndispol
          integer,intent(OUT) :: ret
          integer :: i,j,jj,el,ambi
          logical :: found,allfound,eqname,eqsource
    
          ret  = 0
          ambi = 0
          do i=1,N-1
            if (kind(0).ne.m_kind(0,i)) cycle   ! different no(elements)
            allfound=.true.
            do j=1,kind(0)
              el = kind(j)
              found = .false.
              do jj=1,m_kind(0,i)
                if (el.ne.m_kind(jj,i)) cycle
                found = .true.
                exit
              enddo
              if (.not.found) then
                allfound = .false.
                exit                            ! different elements
              else if (anz(j).ne.m_anz(jj,i)) then
                allfound = .false.
                exit                            ! different stoich.fac.
              endif
            enddo
            if (.not.allfound) cycle
            eqname = (trim(molname)==trim(cmol(i)))
            eqsource = (loop==source(i))
            if (eqname.and.eqsource) then
              print*,"*** double molecule in "//dispolfiles(loop)
              print*,trim(molname)//", "//trim(cmol(i))
              stop
            else if ((.not.eqname).and.eqsource.and.loop==1) then
              if (initchem_info) then
                print*,trim(molname)//", "//trim(cmol(i))// &
                    " different isomere in first source is OK"
              endif
              return  
            else if (eqname.and.(.not.eqsource)) then  
              ret = i
              return
            else
              ambi = i 
            endif
          enddo
          if (ambi>0) then
            if (source(ambi)==loop) then 
              if (initchem_info) then
                print*,trim(molname)//", "//trim(cmol(ambi))// &
                    " different isomere in subsequent source is OK"
              endif  
              ret = 0
              return
            else  
              print*,"*** "//trim(molname)//", "//trim(cmol(ambi))// &
                  " ambiguous names in ..."
              print*,trim(dispolfiles(loop))// &
                  ", "//trim(dispolfiles(source(ambi)))
              print*,"please equalise in both data files."
              stop 
            endif  
          endif  
          end
    
          subroutine copy_molecule_names(nmole,t_cmol)
            use CHEMISTRY,ONLY: cmol
            
            integer, intent(in) :: nmole
            CHARACTER*20, intent(out), dimension(nmole) :: t_cmol
            integer :: ido
            do ido=1,nmole
              t_cmol(ido) = cmol(ido)
            enddo

          end subroutine

          subroutine run_ggchem(nLayers,nnelem,nnmol, elem_out, mol_out)
            use PARAMETERS,ONLY: Tmin,Tmax,pmin,pmax,nHmin,nHmax,useDatabase, &
                                 model_eqcond,model_pconst,Npoints, &
                                 remove_condensates, elements
             use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,catm,el,charge
             use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,muH,amu, &
                                dust_nel,dust_el,dust_nu,dust_nam,dust_mass, &
                                dust_Vol,mass,mel
            use STRUCTURE,ONLY: Npmax,Tgas,press,pelec,dens,nHtot,estruc
             use EXCHANGE,ONLY: nel,nat,nion,nmol,mmol, &
                               H,C,N,O,W,S,Ca,Si,Mg,Al,Fe
             implicit none
             integer,parameter :: qp =16
             integer, intent(in) :: nLayers,nnmol,nnelem
             real*8, intent(out) :: elem_out(nLayers,nNELEM), mol_out(nLayers,nNMOL)
             real*8 :: p,Tg,rhog,rhod,dustV,nHges,nges,mges,kT,pgas
             real*8 :: ff,fold,dmu,dfdmu
             real*8 :: nTEA,pTEA,mu,muold,Jstar,Nstar
             real(kind=qp) :: eps(NELEM),eps00(NELEM)
             real(kind=qp) :: Sat(NDUST),eldust(NDUST),out(NDUST)
             real(kind=qp) :: fac,deps,e_reservoir(NELEM),d_reservoir(NDUST)
             integer :: it,i,j,jj,k,l,NOUT,iW,stindex
             character(len=5000) :: species,NISTspecies,elnames
             character(len=20) :: frmt,name,short_name(NDUST),test1,test2
             character(len=4) :: sup
             character(len=2) :: test3
             character(len=1) :: char
             integer :: verbose=0
             logical :: isOK,hasW,same,TEAinterface=.false.
            
             do i=1,NDUST
              name = dust_nam(i) 
              j=index(name,"[s]")
              short_name(i) = name
              if (j>0) short_name(i)=name(1:j-1)
            enddo
            eps  = eps0
            NOUT = NELM
            if (charge) NOUT=NOUT-1
             e_reservoir = 0.Q0
             d_reservoir = 0.Q0
             eps00 = eps0
             eldust = 0.Q0
            
             muH = 0.0
             do i=1,NELM
                muH = muH + mass(i)*eps(i)
            enddo


             mu = muH
            !  print *,muH
            !  print *,muH/amu
            !  stop
             do i=1,nLayers
                Tg      = Tgas(i)
                p       = press(i) 
                !nHges   = nHtot(i)
                eps0(:) = estruc(i,:)
                nHges = p*mu/(bk*Tg)/muH
                

                do it=1,99
                  nHges = p*mu/(bk*Tg)/muH
                  if (model_eqcond) then
                    call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
                  else
                    eps(:) = eps0(:)
                  endif  
                  !print *, eps
                  !print *,"yes"
                  call GGCHEM(nHges,Tg,eps,.false.,0)
                  kT = bk*Tg
                  nges = nel
                  mges = nel*mel
                  do j=1,NELEM
                    nges = nges + nat(j)
                    mges = mges + nat(j)*mass(j)
                  enddo
                  do j=1,NMOLE
                    nges = nges + nmol(j)
                    mges = mges + nmol(j)*mmol(j)
                  enddo
                  pgas  = nges*bk*Tg
                  ff    = p-pgas
                  if (it==1) then
                    muold = mu
                    mu = nHges/pgas*(bk*Tg)*muH
                    dmu = mu-muold
                  else
                    dfdmu = (ff-fold)/(mu-muold)
                    dmu   = -ff/dfdmu
                    muold = mu
                    if ((dmu>0.0).or.ABS(dmu/mu)<0.7) then
                      mu = muold+dmu
                    else
                      mu = nHges/pgas*(bk*Tg)*muH
                    endif  
                  endif
                  fold = ff
                  !print '("p-it=",i3,"  mu=",2(1pE20.12))',it,mu/amu,dmu/mu
                  if (ABS(dmu/mu)<1.E-10) exit
                enddo
!


            ! --- remove all condensates and put them in reservoir? ---
             if (remove_condensates) then
              fac = 1.Q+0
              do j=1,NDUST
                d_reservoir(j) = d_reservoir(j) + fac*eldust(j)
                do jj=1,dust_nel(j)
                  k = dust_el(j,jj)
                  e_reservoir(k) = e_reservoir(k)  &
                                + fac*dust_nu(j,jj)*eldust(j)
                enddo
              enddo  
              do j=1,NELM
                if (j==el) cycle 
                k = elnum(j)
                !print'(A3,2(1pE18.10))',elnam(k),eps(k)/eps00(k), &
                !              (eps(k)+e_reservoir(k))/eps00(k)
              enddo
              eps0(:) = eps(:) + (1.Q0-fac)*e_reservoir(:)
              estruc(i+1,:) =  eps0(:)  ! for next layer
              eldust(:) = d_reservoir(:)
             endif  
     
             !--- compute supersat ratios and nucleation rates ---
             call SUPERSAT(Tg,nat,nmol,Sat)
             if (hasW) then
              call NUCLEATION('W',Tg,dust_vol(iW),nat(W), &
                            Sat(iW),Jstar,Nstar)
             else
              Jstar = 0
              Nstar = 9.e+30
             endif  
     
             !--- compute dust/gas density ratio ---
             rhog  = nHges*muH
             rhod  = 0.0
             dustV = 0.0
             do jj=1,NDUST
              rhod  = rhod  + nHges*eldust(jj)*dust_mass(jj)
              dustV = dustV + eldust(jj)*dust_Vol(jj)
              out(jj) = LOG10(MIN(1.Q+300,MAX(1.Q-300,Sat(jj))))
              if (ABS(Sat(jj)-1.Q0)<1.E-10) out(jj)=0.Q0
             enddo  


             !print'(i4," Tg[K] =",0pF8.2,"  n<H>[cm-3] =",1pE10.3)', &
             !        i,Tg,nHges

             !write(*,1010) ' Tg=',Tg,' n<H>=',nHges, &
             !                ' p=',pgas/bar,' mu=',mu/amu, &
             !                ' dust/gas=',rhod/rhog
             !elem_out(i,:) = nat

             mol_out(i,:) = real(nmol/nHges,8)

            enddo
            1010 format(A4,0pF8.2,3(a6,1pE9.2),1(a11,1pE9.2))
      end subroutine
end module