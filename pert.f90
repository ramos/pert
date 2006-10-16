
Program Metrica 

  USE NumTypes
  USE Constants
  USE Geometry
  USE Error
  USE Fourier
  USE Bases

  Integer, Parameter :: MaxOrder = 51
  Integer :: Nqs = 2, Norders = 20

  Integer :: Nterm = 20, Terminos(MaxOrder)
  Integer, Allocatable :: Multi(:)
  Complex (kind=DPC), Allocatable :: C(:), Wzero(:)
  Real (kind=DP) :: Dnorm
  Real (kind=DP), Allocatable :: Factor(:)
  Complex (kind=DPC) :: Wfac
  Character (len=100) :: dirbase, FileSave, FileMulti

  Type (Fourier_Serie_2D) :: Chif, Prod, &!For Stage 2
       & Lcont, Aux
  Type (Fourier_Serie_2D), Allocatable :: hf(:), &
       & Pf(:), deltaf(:), Acum2(:), Aparcial(:), &
       & Powhf(:,:) ! For Stage 2

  

  Interface 
     Recursive Function Factorial(N1) Result (Fac)
       USE NumTypes

       Integer, Intent (in) :: N1
       Real (kind=DP) :: Fac
     End Function Factorial
  End Interface

  Write(stderr, *)'Input q, Norders, Nterm, dirbase and C(:)'
  Read(*,*)Nqs, Norders, Nterm, dirbase

  ! Allocate space for all the components...
  Allocate(C(Nqs), Wzero(Nqs))
  Allocate(hf(Norders), Pf(Norders), deltaf(Norders), &
       & Acum2(Norders), Aparcial(Nqs), Powhf(Norders,Norders))
  Allocate(Factor(Norders), Multi(Norders))


  CALL Init_Geometry(qq=Nqs)

  Do I = 1, q
     Read(*,'(2ES33.25)')C(I)
     Write(stderr,'(1I4,2ES20.12)')I, C(I)
  End Do
  Do I = 1, q
     Read(*,'(2ES33.25)')Wzero(I)
     Write(stderr,'(1I4,2ES20.12)')I, Wzero(I)
  End Do

  Dnorm = Sqrt(Sum(Abs(C(:))**2))
  C = C/Dnorm

  ! Set values to Terminos(51)
  Data Terminos /0, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101,&
       & 135, 176, 231, 297, 385, 490, 627, 792, 1002, 1255, 1575, &
       & 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, &
       & 12310, 14883, 17977, 21637, 26015, 31185, 37338, 44583, &
       & 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525, &
       & 204226/
  

  Write(0,*)
  Write(0,*)'Stage 1: Calculating h, Factor:'
  Write(0,*)'==============================='

  ! Init the fourier Series
  CALL Init_Serie(Chif, Nterm)
  CALL Init_Serie(Prod, Nterm)
  
  Do I = 1, q
     CALL Init_Serie(Aparcial(I), Nterm)
  End Do
    
  Do I = 1, Norders
     CALL Init_Serie(hf(I), Nterm)
     CALL Init_Serie(Deltaf(I), Nterm)
     CALL Init_Serie(Acum2(I), Nterm)
     CALL Init_Serie(Pf(I), Nterm)
     Do K = 1, Norders
        CALL Init_Serie(Powhf(I,K), Nterm)
     End Do
  End Do

  ! Set the displacements and
  ! calculate the initial function.
  Do N1 = -Nterm, Nterm
     Do N2 = -Nterm, Nterm
        Chif%Coef(N1, N2) = (0.0_DP, 0.0_DP)
        Do I = 1, q
           Do J = 1, q
              Chif%Coef(N1,N2) = Chif%Coef(N1,N2) + &
                   & L(i, j, N1, N2)*Conjg(C(i))*C(j)
           End Do
        End Do
     End Do
  End Do
  
  I = 1
  ! Set "maually" the first order.
  Factor(1) = Real(1.0_DP/Chif%Coef(0,0), kind=DP)
  
  Deltaf(1) = -0.5_DP * Factor(1) * Chif
  Deltaf(1)%Coef(0,0) = (0.0_DP, 0.0_DP)

  Pf(1)%Coef = (0.0_DP, 0.0_DP)
  Pf(1)%Coef(0,0) = Cmplx(Factor(1), kind=DPC)

  
  Do N1 = -Nterm, Nterm
     Do N2 = -Nterm, Nterm
        If ( (N1 == 0) .and. (N2 == 0) ) Then
           hf(1)%Coef(0,0) = (0.0_DP, 0.0_DP)
        Else
           hf(1)%Coef(N1,N2) = - Deltaf(1)%Coef(N1,N2) / xi(N1,N2)
        End If
     End Do
  End Do
  
  Write(FileSave, '(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'factor:O=',I&
       &,':Nterm=',Nterm ,':flux=',q,'.dat'
  FileSave = Trim(Trim(dirbase) // '/' // FileSave)
  Open (Unit=99, File = FileSave)
  Write(99,'(1ES33.25)')Factor(1)
  Close(99)
  
  Write(FileSave,'(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'deltaf:O=',I&
       &,':Nterm=',Nterm ,':flux=',q,'.dat'
  FileSave = Trim(Trim(dirbase) // '/' // FileSave)
  CALL Save_Serie(Deltaf(I), Trim(FileSave))
  
  Write(FileSave,'(1A5,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'hf:O=',I&
       &,':Nterm=',Nterm ,':flux=',q,'.dat'
  FileSave = Trim(Trim(dirbase) // '/' // FileSave)
  CALL Save_Serie(hf(I), Trim(FileSave))
  
  Powhf(1,1) = hf(1)
  Do I2 = 2, Int((Norders-1))
     Powhf(1,I2) = hf(1) * Powhf(1,I2-1)
  End Do
  
  
  ! Now start the iteration until Norders is reached.
  Do I = 2, Norders
     ! This loop must be done for each
     ! order. I is the order.
     Write(FileMulti, '(1A15,1I2.2,1A4)')'combinatoria/O=', I-1, '.dat'     
     Open (Unit = 34, File = Trim(FileMulti), ACTION="READ")     
     
     Acum2(I)%Coef = (0.0_DP, 0.0_DP)
     Do J = 1, Terminos(I)
        Read(34,*)(Multi(K), K=1, Norders-1)
        Ntot = Sum(Multi)
        
        CALL Unit(Prod,Nterm)
        Do K = 1, I-1
           If (Multi(K) /= 0) Then
              Prod = (1.0_DP/Factorial(Multi(K))) * Prod * &
                   & Powhf(K,Multi(K))
           End If
        End Do
        Prod = (-2.0_DP) ** Ntot * Prod
        
        Acum2(I) = Acum2(I) + Prod
     End Do
     close(34)
     
     Prod%Coef = (0.0_DP, 0.0_DP)
     Do J = 1, I-1
        Prod = Prod + Factor(I-J) * Acum2(J+1)
     End Do
     
     Pf(I) = Prod
     Prod = Prod * Chif
     
     Factor(I) = -Real(Prod%Coef(0,0),kind=DP) / Real(Chif%Coef(0,0),kind=DP)
     Deltaf(I) = -0.5_DP * ( (Factor(I)*Chif) + Prod )

     Pf(I)%Coef(0,0) = Pf(I)%Coef(0,0) + Factor(I)
     Write(FileSave,'(1A5,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'Pf:O=',I&
          &,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     CALL Save_Serie(Pf(I), Trim(FileSave))

     
     Do N1 = -Nterm, Nterm
        Do N2 = -Nterm, Nterm
           If ( (N1 == 0) .and. (N2 == 0) ) Then
              hf(I)%Coef(0,0) = (0.0_DP, 0.0_DP)
           Else
              hf(I)%Coef(N1,N2) = hf(I-1)%Coef(N1,N2) - &
                   & Deltaf(I)%Coef(N1,N2) / xi(N1,N2)
           End If
        End Do
     End Do
     
     Write(FileSave, '(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'factor:O=',I&
          &,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     Open (Unit=99, File = FileSave)
     Write(99,'(1ES33.25)')Factor(I)
     Close(99)
     
     Write(FileSave,'(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'deltaf:O=',I&
          &,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     CALL Save_Serie(Deltaf(I), Trim(FileSave))
     
     Write(FileSave,'(1A5,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'hf:O=',I&
          &,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     CALL Save_Serie(hf(I), Trim(FileSave))
     
     Powhf(I,1) = hf(I)
     Do Npow = 2, Int((Norders-1)/I)
        Powhf(I, Npow) = hf(I) * Powhf(I,Npow-1)
     End Do
     
     Write(stderr, *)I, Deltaf(I)%Coef(0,0)
  End Do


  Stop
End Program Metrica

! ****************************
! *
Recursive Function Factorial(N1) Result (Fac)
! *
! ****************************
! * Returns the factorial of a 
! * integer number.
! ****************************
  
  USE NumTypes

  Integer, Intent (in) :: N1
  Real (kind=DP) :: Fac
  
  
  If (N1 == 0) Then 
     Fac = 1.0_DP
  Else
     Fac = Real(N1,kind=DP) * Factorial(N1 - 1)
  End If
  
  Return
End Function Factorial
