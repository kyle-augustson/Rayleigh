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

Include 'mkl_sparse_qr.f90'

Program ULTRA_Unit_Tests
  Use Ultraspherical_Polynomials
  Use mkl_sparse_qr
  Implicit None

  Call Test_Derivatives()
  
Contains

  Subroutine Test_Derivatives()
    Implicit None
    Integer, Parameter :: N=64
    Integer :: i, lambda, stat
    Real*8 :: pi, x(N), ifx(N,1), fx(N,1), coefs(N,1), tcoefs(N,1), sol(N,1)
    type(sparse_matrix_t) :: Dl, C0, C1, C2
    type(matrix_descr) :: descr

    ! There could be some tuning of the descriptor to select the appropriate upper
    ! triangular matrix parameters for faster operations, when appropriate (Derivative and Conversion)
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    
    pi = 3.141592653589793238
    
    Do i=0,N-1
       x(i+1) = -cos(pi*(i+0.5)/N)
    End Do

    fx(:,1) = x**8
    Call Forward_Chebyshev_Transform(x,fx,coefs)
    Call Backward_Chebyshev_Transform(x,coefs,ifx)
    Print*, 'L1 Chebyshev Transform Error'
    Print*, maxval(abs(fx(:,1)-ifx(:,1)))

    Call Build_Conversion_Matrix(N,0,C0,stat)
    
    Do lambda=1,2
       Call Build_Derivative_Matrix(N,lambda,Dl,stat)
       stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,1d0,Dl,descr,coefs,0d0,tcoefs)
       !Print*, lambda, tcoefs
       !stat = mkl_sparse_d_qr(SPARSE_OPERATION_NON_TRANSPOSE,C0,descr,&
       !     SPARSE_LAYOUT_COLUMN_MAJOR,1,sol,N,tcoefs,N)
       ! This solve seems like garbage may be just go dense for now and use dgetrf ?
       ! Or just use basis recombination from the get go and banded solve?
       Call dgetrf(N,N,forward_transform_matrix,N,pivots,info))
       pi = 1
       Do i=0,lambda-1
          pi = pi*(8-i)
       End Do
       fx(:,1) = pi*x**(8-lambda)
       Call Forward_Chebyshev_Transform(x,fx,tcoefs)
       Print*, 'Derivative order lambda'
       Print*, lambda
       Print*, 'L1 coefficient error'
       Print*, maxval(abs(sol-tcoefs))
       Call Build_Conversion_Matrix(N,lambda,C1,stat)
       stat = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE,descr,C1, &
            SPARSE_OPERATION_NON_TRANSPOSE,descr,C0,SPARSE_STAGE_FULL_MULT,C2)
       stat = mkl_sparse_copy(C2,descr,C0)
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
  
End Program ULTRA_Unit_Tests
  

  
