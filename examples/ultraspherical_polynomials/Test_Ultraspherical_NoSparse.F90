!
!  A library in Python, C, and Fortran for Ultraspherical polynomials
!  and sparse spectral methods for PDEs (ULTRA).
!  
!  Copyright (C) 2023 Kyle Augustson.
!
!  This library is distributed along side RAYLEIGH
!  (Copyright (C) 2018 the authors of RAYLEIGH).
!
!  ULTRA and RAYLEIGH are free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  ULTRA and RAYLEIGH are distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with ULTRA and RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!

!  This program tests the Fortran subroutines of the Ultraspherical polynomial library.

Program ULTRA_Unit_Tests_NoSparse
  Use Ultraspherical_Polynomials_NoSparse
  Implicit None

  Call Test_Derivatives()
  
Contains

  Subroutine Test_Derivatives()
    Implicit None
    Integer, Parameter :: N=64
    Integer :: i, j, k, lambda, info, pivots(N)
    Real*8 :: pi, x(N), ifx(N,1), fx(N,1), coefs(N,1), tcoefs(N,1), sol(N,1), exact(N,8)
    Real*8, Dimension(N,N) :: Dl, C0, C1, C2, LU

    pi = 3.141592653589793238d0
    
    Do i=0,N-1
       x(i+1) = -cos(pi*(Dble(i)+0.5d0)/Dble(N))
    End Do

    !T_8(x)
    fx(:,1) = 128d0*x**8 - 256d0*x**6 + 160d0*x**4 - 32d0*x**2 + 1d0

    exact(:,1) = 8d0*128d0*x**7 - 6d0*256d0*x**5 + 4d0*160d0*x**3 - 2d0*32d0*x
    exact(:,2) = 7d0*8d0*128d0*x**6 - 5d0*6d0*256d0*x**4 + 3d0*4d0*160d0*x**2 - 2d0*32d0
    exact(:,3) = 6d0*7d0*8d0*128d0*x**5 - 4d0*5d0*6d0*256d0*x**3 + 2d0*3d0*4d0*160d0*x
    exact(:,4) = 5d0*6d0*7d0*8d0*128d0*x**4 - 3d0*4d0*5d0*6d0*256d0*x**2 + 2d0*3d0*4d0*160d0
    exact(:,5) = 4d0*5d0*6d0*7d0*8d0*128d0*x**3 - 2d0*3d0*4d0*5d0*6d0*256d0*x
    exact(:,6) = 3d0*4d0*5d0*6d0*7d0*8d0*128d0*x**2 - 2d0*3d0*4d0*5d0*6d0*256d0
    exact(:,7) = 2d0*3d0*4d0*5d0*6d0*7d0*8d0*128d0*x
    exact(:,8) = 2d0*3d0*4d0*5d0*6d0*7d0*8d0*128d0
    
    Call Forward_Chebyshev_Transform(x,fx,coefs)
    Print*, 'T_n Coefs for f'
    Print*, coefs
    tcoefs = 0d0
    tcoefs(9,1) = 1d0
    Print*, 'Linf Forward error'
    Print*, maxval(abs(coefs-tcoefs))
    
    Call Backward_Chebyshev_Transform(x,tcoefs,ifx)
    Print*, 'Linf Backward error'
    Print*, maxval(abs(fx(:,1)-ifx(:,1)))/maxval(abs(fx))

    Call Backward_Chebyshev_Transform(x,coefs,ifx)
    Print*, 'Linf Forward-Backward error'
    Print*, maxval(abs(fx(:,1)-ifx(:,1)))/maxval(abs(fx))
    
    Call Build_Conversion_Matrix(N,0,C0)
    
    Do lambda=1,8
       Call Build_Derivative_Matrix(N,lambda,Dl)
       tcoefs = 0d0
       Do i=1,N
          Do j=1,N
             tcoefs(i,1) = tcoefs(i,1) + Dl(i,j)*coefs(j,1)
          End Do
       End Do

       fx(:,1) = exact(:,lambda)
       Call Forward_Chebyshev_Transform(x,fx,sol)
       ! Convert to C^\lambda basis
       ifx = 0d0
       Do i=1,N
          Do j=1,N
             ifx(i,1) = ifx(i,1) + C0(i,j)*sol(j,1)
          End Do
       End Do
       Print*, 'Derivative order lambda'
       Print*, lambda
       Print*, 'Linf coefficient error'
       Print*, maxval(abs(ifx-tcoefs))/maxval(abs(ifx))
       Call Build_Conversion_Matrix(N,lambda,C1)
       C2 = 0
       Do j=1,N
          Do i=1,N
             Do k=1,N
                C2(i,j) = C2(i,j) + C1(i,k)*C0(k,j)
             End Do
          End Do
       End Do
       C0 = C2
    End Do
  End Subroutine Test_Derivatives

  !Subroutine Test_Multiplication()
  !  Implicit None
  !  Integer, Parameter :: N=64
  !  Integer :: i, lambda, stat
  !  Real*8 :: x(N), fx(N), coefs(N)
  !  type(sparse_matrix_t) :: Dl, C(8)
  !
  !  Do i=0,N-1
  !  End Do
  !  Call Build_Conversion_Matrix(N,0,C(0),stat)
  !  Do lambda=1,8
  !     Call Build_Conversion_Matrix(N,0,C(lambda),stat)
  !     Call Build_Derivative_Matrix(N,lambda,matrix,stat)
  !  End Do
  !End Subroutine Test_Multiplication

  !Subroutine Test_Boundaries()
  !End Subroutine Test_Boundaries

  !Subroutine Test_Problem_IVP()
  !End Subroutine Test_Problem_IVP

  !Subroutine Test_Problem_BVP()
  !End Subroutine Test_Problem_BVP
  
End Program ULTRA_Unit_Tests_NoSparse
  

  
