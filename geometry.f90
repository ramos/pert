! ******************************************
! *
MODULE Geometry
! *
! ******************************************
! *
! * Data with the geometry of the Torus.
! *
! ******************************************
 
  USE NumTypes
  USE Constants

  Integer :: q
  Real (kind=DP) :: l1, l2, taui, epsilon
  Complex (kind=DPC) :: tau


CONTAINS

! ******************************************
! *
  Subroutine Init_Geometry(qq, ll1, ll2)
! *
! ******************************************
    
    Real (kind=DP), Optional :: ll1, ll2
    Integer, Optional :: qq

    If (Present(qq)) Then
       q = qq
    Else 
       q = 1   
    End If

    If ( Present(ll1).and.Present(ll2) ) Then
       l1 = ll1
       l2 = ll2
    Else 
       taui = 1.0_DP
       l2 = Sqrt(taui * 4.0_DP * q * PI_DP/(1.0_DP-0.5_DP))
       l1 = Sqrt(4.0_DP * q * PI_DP/(taui* (1.0_DP-0.5_DP)))
    End If

    taui = l2/l1
    tau = UNITIMAG_DPC*taui
    epsilon = 1.0_DP - 4.0_DP*q*PI_DP/(l1*l2)


    Return
  End Subroutine Init_Geometry
    
! ******************************************
! *
  Real (kind=DP) Function Critical_Area()
! *
! ******************************************
    
    Critical_Area = 4.0_DP * PI_DP * q

    Return
  End Function Critical_Area

End MODULE Geometry
