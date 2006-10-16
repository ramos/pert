
! ******************************************
! *
MODULE Bases
! *
! ******************************************
! *
! * Contains the definition of the basis of
! * the moduli space of solutions.
! *
! ******************************************
 
  USE NumTypes
  USE Constants
  USE Specialfunc
  USE Error
  
  USE Geometry


  Integer, Parameter :: NF = 20

  Private NF

CONTAINS

! ******************************************
! *
  Complex (kind=DPC) Function psi(i, z, Tol)
! *
! ******************************************
    
    Integer, Intent (in) :: i
    Complex (kind=DPC) :: z
    Real (kind=DP), Optional :: Tol

    Real (kind=DP) :: Err
    
    If (Present(Tol)) Then
       psi = Sqrt(Sqrt(2.0_DP*q*taui)) * &
            & ThetaChar(Real(i,kind=DP)/Real(q,kind=DP), 0.0_DP, &
            & TWOPI_DP*z*q/l1, q*tau, Tol)
    Else       
       psi = Sqrt(Sqrt(2.0_DP*q*taui)) * &
            & ThetaChar(Real(i,kind=DP)/Real(q,kind=DP), 0.0_DP, &
            & TWOPI_DP*z*q/l1, q*tau)
    End If


    Return
  End Function psi

! ******************************************
! *
  Function Cbase(w) Result (c)
! *
! ******************************************
    
    Complex (kind=DPC), Intent (in) :: w(:)
    Complex (kind=DPC) :: c(Size(w,1))

    Real (kind=DP) :: Dnorm

    If (q /= Size(w,1)) Then
       CALL Abort("Cbase", "Topological Charge (q) and number of zeros must coincide.")
    End If


    If (q == 2) Then
       c(1) = psi(1, w(1))
       c(2) = -psi(0, w(1))
    Else 
       c = UNITIMAG_DPC
    End If

    ! Normalise the c's
    Dnorm = 0.0_DP
    Do I = 1, q
       Dnorm = Dnorm + Abs(c(I))**2
    End Do
    Dnorm = Sqrt(Dnorm)

    Do I = 1, q
       c(I) = c(I)/Dnorm
    End Do

    Return
  End Function Cbase

! **************************************
! *
  Function xi(n1, n2)
! *
! **************************************
    
    Integer, Intent (in) :: n1, n2
    Real (kind=DP) :: xi
    
    xi = PI_DP * taui / q * &
         & ( Real(n1,kind=DP)**2 + Real(n2,kind=DP)**2/taui**2 )
    
    Return
  End Function xi

! **************************************
! *
  Function L(i, j, n1, n2)
! *
! **************************************

    Integer, Intent (in) :: i, j, n1, n2
    Complex (kind=DPC) :: L
    
!    If ((i-1) == (Mod(j+n1-1,q))) Then
    If ( Mod( (i-1) - (Mod(j+n1-1,q)), q ) == 0) Then
       L = exp(-TWOPI_IMAG_DPC/q * ((j-1)*n2 + &
            & Real(n1*n2,kind=DP)/2.0_DP) - xi(n1, n2)/2.0_DP ) 
    Else
       L = (0.0_DP, 0.0_DP)
    End If
    
    Return
  End Function L


! **************************************
! *
  Function ZeroRes(n, Func, l1, l2, TTol)
! *
! **************************************

    USE Integration

    Real (kind=DP), Intent (in) :: l1, l2
    Real (kind=DP), Intent (in), Optional :: TTol
    Integer, Intent (in) :: n
    Complex (kind=DPC) :: ZeroRes

    Real (kind=DP) :: Tol = 1.E-2_DP

    Interface
       Function Func(z)
         USE NumTypes
         
         Complex (kind=DPC), Intent (in) :: z
         Complex (kind=DPC) :: Func
       End Function Func
    End Interface

    ZeroRes = (0.0_DP, 0.0_DP)


    Return
  End Function ZeroRes
  

! **************************************
! *
  Function Fourier3(z)
! *
! **************************************
    
    
    Complex (kind=DPC), Intent(in) :: Z
    Complex (kind=DPC) :: Fourier3, Term
    
    Real (kind=DP) :: X, Y
    
    
    X = Real(z, kind=DP)
    Y = Real(Aimag(z), kind=DP)
    Fourier3 = (0.D0, 0.D0)
    Do I = -NF, NF
       Do J = -NF, NF
          Term =  exp(-PI_DP*q/2.0_DP*taui * &
               & (Real(I**2, kind=DP)+(Real(J,kind=DP)/taui)**2)) * &
               & (-1)**(I*J) * &
               & exp(TWOPI_IMAG_DPC*q*X*I/l1) * &
               & exp(TWOPI_IMAG_DPC*q*Y/l2*J/taui)
          
          Fourier3 = Fourier3 + Term
       End Do
    End Do
    
    Return
  End Function Fourier3

! **************************************
! *
  Function Fourier1(z)
! *
! **************************************
    
    
    Complex (kind=DPC), Intent(in) :: Z
    Complex (kind=DPC) :: Fourier1, Term
    
    Real (kind=DP) :: X, Y
    
    
    X = Real(z, kind=DP)
    Y = Real(Aimag(z), kind=DP)
    Fourier1 = (0.D0, 0.D0)
    Do I = -NF, NF
       Do J = -NF, NF
          Term =  exp(-PI_DP*q/2.0_DP*taui * &
               & (Real(I**2, kind=DP)+(Real(J,kind=DP)/taui)**2)) * &
               & (-1)**(I*J) * (-1)**(I+q*J) * &
               & exp(TWOPI_IMAG_DPC*q*X*I/l1) * &
               & exp(TWOPI_IMAG_DPC*q*Y/l2*J/taui)
          
          Fourier1 = Fourier1 + Term
       End Do
    End Do
    
    Return
  End Function Fourier1

! **************************************
! *
  Function Fourier4(z)
! *
! **************************************
    
    
    Complex (kind=DPC), Intent(in) :: Z
    Complex (kind=DPC) :: Fourier4, Term
    
    Real (kind=DP) :: X, Y
    
    
    X = Real(z, kind=DP)
    Y = Real(Aimag(z), kind=DP)
    Fourier4 = (0.D0, 0.D0)
    Do I = -NF, NF
       Do J = -NF, NF
          Term =  exp(-PI_DP*q/2.0_DP*taui * &
               & (Real(I**2, kind=DP)+(Real(J,kind=DP)/taui)**2)) * &
               & (-1)**(I*J) * (-1)**I * &
               & exp(TWOPI_IMAG_DPC*q*X*I/l1) * &
               & exp(TWOPI_IMAG_DPC*q*Y/l2*J/taui)
          
          Fourier4 = Fourier4 + Term
       End Do
    End Do
    
    Return
  End Function Fourier4

! **************************************
! *
  Subroutine Fourier1234(z, b, Fourier1, Fourier2, Fourier3, Fourier4)
! *
! **************************************
    
    
    Complex (kind=DPC), Intent(in)  :: Z
    Real (kind=DP) :: b
    Real (kind=DP), Intent(out) :: Fourier1, Fourier2, &
         & Fourier3, Fourier4
    
    Complex (kind=DPC) :: Term, F1, F2, F3, F4, Factor
    Real (kind=DP) :: X, Y
    
    
    X = Real(z, kind=DP)
    Y = Real(Aimag(z), kind=DP)
    F1 = 0.D0
    F2 = 0.D0
    F3 = 0.D0
    F4 = 0.D0
    Do I = -NF, NF
       Do J = -NF, NF
!          Term =  exp(-PI_DP*q/2.0_DP*taui * &
!               & (Real(I**2, kind=DP)+(q*Real(J,kind=DP)/taui)**2)) * &
!               & (-1)**(I*J) * &
!               & exp(TWOPI_IMAG_DPC*q*X*I/l1) * &
!               & exp(TWOPI_IMAG_DPC*q*Y/l2*J/taui)
          
!          F1 = F1 + Term * (-1)**(I+q*J)
!          F2 = F2 + Term * (-1)**(q*J)  
!          F3 = F3 + Term  
!          F4 = F4 + Term * (-1)**I
          Term = exp(-HALFPI_DP*b * &
               & (Real(I,kind=DP)**2 + (Real(J,kind=DP)/b)**2) ) * &
               & exp(2.0_DP*UNITIMAG_DPC*(Real(I,kind=DP)*X + Real(J&
               &,kind=DP)*Y/b))

          F1 = F1 + Term * (-1.0_DP)**(I*J + I+J)
          F2 = F2 + Term * (-1.0_DP)**(I*J + J)
          F3 = F3 + Term * (-1.0_DP)**(I*J)
          F4 = F4 + Term * (-1.0_DP)**(I*J + I)


       End Do
    End Do

    Factor = exp(2.0_DP*Y**2/(PI_DP*b))/Sqrt(2.0_DP*b)
    Fourier1 = Real(F1,kind=DP) * Factor 
    Fourier2 = Real(F2,kind=DP) * Factor 
    Fourier3 = Real(F3,kind=DP) * Factor 
    Fourier4 = Real(F4,kind=DP) * Factor 
    
    Return
  End Subroutine Fourier1234


End MODULE Bases
