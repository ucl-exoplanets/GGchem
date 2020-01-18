module fort_ggchem
  implicit none



  contains
      SUBROUTINE INIT_LEAN()
        use DUST_DATA,ONLY: NELEM,eps0,mass,muH,elnam,amu
        use EXCHANGE,ONLY: H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl, & 
                          Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge, & 
                           As,Se,Br,Kr,Rb,Sr,Y,Zr,W 
      implicit none
      integer,parameter :: qp=selected_real_kind(33,4931)
      integer :: i,j,nr
      real(kind=qp) :: abund(74,4),eps(NELEM),epsH,mfrac(NELEM)
      real(kind=qp) :: m,val,addH2O
      character(len=2) :: el
      character(len=20) :: elname
      character(len=10) :: source(4)
      character(len=200) :: line
      logical :: found

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


      SUBROUTINE INIT_TAUREX_CHEMISTRY(nelem_t, cel, ndispol, dispolfiles, do_charge)


        use CHEMISTRY,ONLY: NMOLdim,NMOLE,NELM,catm,cmol,el, &
           source,fit,natom,a,error,i_nasa, &
            m_kind,m_anz,elnum,elion,charge, &
            el,H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc, &
            Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,W
        use DUST_DATA,ONLY: mass,mel,amu
        use EXCHANGE,ONLY: nmol,mmol

        implicit none

        character*1024,intent(in) :: dispolfiles(ndispol)
        character*2, intent(in) :: cel(nelem_t)
        character*2 :: elnam
        integer, intent(in) :: nelem_t, ndispol
        logical, intent(in) :: do_charge

        integer :: loop,i,ii,j,iel,e,smax,ret

        NELM = nelem_t

        charge = do_charge

        do i=1,nelem_t
                elnam = cel(i)
                catm(i) = elnam
                print*,'element '//elnam
                if     (elnam=='el') then; el=i ; charge=.true.
                elseif (elnam=='H')  then;  H=i ; elnum(i)=1 
                elseif (elnam=='He') then; He=i ; elnum(i)=2 
                elseif (elnam=='Li') then; Li=i ; elnum(i)=3
                elseif (elnam=='Be') then; Be=i ; elnum(i)=4 
                elseif (elnam=='B')  then;  B=i ; elnum(i)=5 
                elseif (elnam=='C')  then;  C=i ; elnum(i)=6 
                elseif (elnam=='N')  then;  N=i ; elnum(i)=7 
                elseif (elnam=='O')  then;  O=i ; elnum(i)=8 
                elseif (elnam=='F')  then;  F=i ; elnum(i)=9 
                elseif (elnam=='Ne') then; Ne=i ; elnum(i)=10
                elseif (elnam=='Na') then; Na=i ; elnum(i)=11 
                elseif (elnam=='Mg') then; Mg=i ; elnum(i)=12
                elseif (elnam=='Al') then; Al=i ; elnum(i)=13
                elseif (elnam=='Si') then; Si=i ; elnum(i)=14
                elseif (elnam=='P')  then;  P=i ; elnum(i)=15 
                elseif (elnam=='S')  then;  S=i ; elnum(i)=16
                elseif (elnam=='Cl') then; Cl=i ; elnum(i)=17
                elseif (elnam=='Ar') then; Ar=i ; elnum(i)=18
                elseif (elnam=='K')  then;  K=i ; elnum(i)=19
                elseif (elnam=='Ca') then; Ca=i ; elnum(i)=20
                elseif (elnam=='Sc') then; Sc=i ; elnum(i)=21
                elseif (elnam=='Ti') then; Ti=i ; elnum(i)=22
                elseif (elnam=='V')  then;  V=i ; elnum(i)=23
                elseif (elnam=='Cr') then; Cr=i ; elnum(i)=24
                elseif (elnam=='Mn') then; Mn=i ; elnum(i)=25
                elseif (elnam=='Fe') then; Fe=i ; elnum(i)=26
                elseif (elnam=='Co') then; Co=i ; elnum(i)=27
                elseif (elnam=='Ni') then; Ni=i ; elnum(i)=28
                elseif (elnam=='Cu') then; Cu=i ; elnum(i)=29
                elseif (elnam=='Zn') then; Zn=i ; elnum(i)=30
                elseif (elnam=='Ga') then; Ga=i ; elnum(i)=31
                elseif (elnam=='Ge') then; Ge=i ; elnum(i)=32 
                elseif (elnam=='As') then; As=i ; elnum(i)=33 
                elseif (elnam=='Se') then; Se=i ; elnum(i)=34 
                elseif (elnam=='Br') then; Br=i ; elnum(i)=35 
                elseif (elnam=='Kr') then; Kr=i ; elnum(i)=36 
                elseif (elnam=='Rb') then; Rb=i ; elnum(i)=37 
                elseif (elnam=='Sr') then; Sr=i ; elnum(i)=38 
                elseif (elnam=='Y')  then;  Y=i ; elnum(i)=39 
                elseif (elnam=='Zr') then; Zr=i ; elnum(i)=40
                elseif (elnam=='W')  then;  W=i ; elnum(i)=41
                else
                stop "*** unknown element "
                endif
                !    if (initchem_info) then
                print*,'element '//elnam,elnum(NELM)
                !    endif
            !enddo
          enddo

          end subroutine
end module