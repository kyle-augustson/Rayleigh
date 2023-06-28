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
    Integer, Parameter :: N=512
    Integer :: i, j, k, lambda, info, lwork_qr, pivots(N)
    Real*8 :: pi, x(N), ifx(N,1), fx(N,1), coefs(N,1), tcoefs(N,1), sol(N,1), Pre(N), tau_qr(N)
    Real*8, Allocatable :: work_qr(:)
    Real*8, Dimension(N,N) :: Dl, C0, C1, C2, LU

    pi = 3.141592653589793238d0
    
    Do i=1,N
       x(i) = cos(pi*(Dble(i)-0.5d0)/Dble(N))
    End Do

    Do lambda=1,8
       Call Build_Derivative_Matrix(N,lambda,Dl)

       ! Preconditioner
       Do i=1,N
          Pre(i) = 1d0/Dble(i+lambda-1)
       End Do
       Pre = Pre/(2d0**(lambda-1)*gamma(Dble(lambda)))

       tcoefs = 0d0
       Do i=1,N
          tcoefs(:,1) = tcoefs(:,1) + Dl(:,i)*coefs(i,1)
       End Do

       ! Solve for T_n coefs (invert conversion)
       sol = 0d0
       LU = 0d0
       Do i=1,N
          sol(i,1) = Pre(i)*tcoefs(i,1)
          LU(:,i) = Pre(:)*C0(:,i)
       End Do
       C1 = LU
       Call dgetrf(N,N,LU,N,pivots,info)
       C2 = 0d0
       Do i=1,N
          C2(i,i) = 1d0
       End Do
       Print*, 'Derivative order lambda'
       Print*, lambda
       Call dgetrs('N',N,N,LU,N,pivots,C1,N,info)
       Print*, 'Linf LU solve C0^(-1) C0 ~ I'
       Print*, maxval(abs(C1-C2))
       Call dgetrs('N',N,1,LU,N,pivots,sol,N,info)
       Print*, 'Linf coefficient error using LU'
       Print*, maxval(abs(exact_coef(:,lambda)-sol(:,1)))/maxval(abs(exact_coef(:,lambda)))
       Call Backward_Chebyshev_Transform(x,sol,ifx)
       Print*, 'Linf spatial error using LU'
       Print*, maxval(abs(exact(:,lambda)-ifx(:,1)))/maxval(abs(exact(:,lambda)))

       sol = 0d0
       LU = 0d0
       Do i=1,N
          sol(i,1) = Pre(i)*tcoefs(i,1)
          LU(:,i) = Pre(:)*C0(:,i)
       End Do
       Allocate(work_qr(N))
       Call dgeqrf(N,N,LU,N,tau_qr,work_qr,-1,info)
       lwork_qr = int(work_qr(1))
       Deallocate(work_qr)
       Allocate(work_qr(lwork_qr))
       sol = 0d0
       LU = 0d0
       Do i=1,N
          sol(i,1) = Pre(i)*tcoefs(i,1)
          LU(:,i) = Pre(:)*C0(:,i)
       End Do
       Call dgeqrf(N,N,LU,N,tau_qr,work_qr,lwork_qr,info)
       Call dormqr('L','T',N,1,N,LU,N,tau_qr,sol,N,work_qr,lwork_qr,info)
       Call dtrsm('L','U','N','N',N,1,1d0,LU,N,sol,N)
       Print*, 'Linf coefficient error using QR'
       Print*, maxval(abs(exact_coef(:,lambda)-sol(:,1)))/maxval(abs(exact_coef(:,lambda)))
       Call Backward_Chebyshev_Transform(x,sol,ifx)
       Print*, 'Linf spatial error using QR'
       Print*, maxval(abs(exact(:,lambda)-ifx(:,1)))/maxval(abs(exact(:,lambda)))
       Deallocate(work_qr)
       Call Build_Conversion_Matrix(N,lambda,C1)
       C2 = 0
       Do j=1,N
          Do k=1,N
             C2(:,j) = C2(:,j) + C1(:,k)*C0(k,j)
          End Do
       End Do
       C0 = C2
    End Do
  End Subroutine Test_Derivatives
  
End Program ULTRA_Unit_Tests_NoSparse
  

  
