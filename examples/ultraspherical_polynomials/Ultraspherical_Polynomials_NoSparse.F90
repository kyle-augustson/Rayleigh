!
!  A library in Python, C, and Fortran for Ultraspherical polynomials
!  and sparse spectral methods for PDEs (ULTRA).
!  
!  Copyright (C) 2023 Kyle Augustson.
!
!  This library is distributed along side RAYLEIGH.
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

! This module defines the matrices necessary to construct sparse
! spectral operators using Ultraspherical polynomials.
! This file thus includes routines to build the following:
! 1. The derivative matrices of arbitrary order.
! 2. The conversion matrices lifting Chebyshev coefficients defined
!    over the Chebyshev polynomials of the first kind to various
!    orders of Ultraspherical polynomials 
! 3. Multiplication matrices for non-constant coefficients that
!    link nccs to a field.
! 4. Forward and backward Chebyshev transform for a nearly arbitrary grid
!    (depends on Lapack to compute LU decomposition of the backward transform,
!      which can fail depending on the structure of the grid)

Module Ultraspherical_Polynomials_NoSparse
  Implicit None

  ! Alloctable arrays for the transforms, only allocated upon the first calls
  ! to the forward or backward transforms.
  Integer, Allocatable :: pivots(:)
  Real*8, Allocatable :: forward_transform_matrix(:,:)
  Real*8, Allocatable :: backward_transform_matrix(:,:)
  
Contains

  ! Building a transform matrix allows for arbitrary grids.
  ! Do the forward transform by doing an LU decomposition of the
  ! backward transform and then solving for the coefs.
  Subroutine Build_Forward_Transform(x)
    Implicit None
    Real*8, Intent(In) :: x(:)
    Integer :: N, info
    If (.not. Allocated(forward_transform_matrix)) Then
       N = size(x)
       Allocate(forward_transform_matrix(N,N),pivots(N))
       Call Build_Backward_Transform(x)
       forward_transform_matrix = backward_transform_matrix
       Call dgetrf(N,N,forward_transform_matrix,N,pivots,info)
       ! Trap error
       If (info<0) Then
          Print*, 'Illegal argument to dgetrf at ', info
          Print*, 'Check input grid values.'
       Else If(info>0) Then
          Print*, 'Factorization for forward Chebyshev transform is singular due to element at ', info
          Print*, 'Adjust the input grid values.'
       End If
    End If
  End Subroutine Build_Forward_Transform


  ! Building a transform matrix allows for arbitrary grids at the expense
  ! of more multiplication operators compared to the DCT.
  ! Explicitly define the backward transform using the properties of the
  ! Chebyshev polynomials of the first kind.
  Subroutine Build_Backward_Transform(x)
    Implicit None
    Real*8, Intent(In) :: x(:) ! The grid points.
    Real*8, Allocatable :: tmp(:), grid(:)
    Real*8 :: x0, x1
    Integer :: N, i, j
    ! On first entry allocate and define the matrix.
    If (.not. Allocated(backward_transform_matrix)) Then
       N = size(x)
       Allocate(backward_transform_matrix(N,N),tmp(N),grid(N))

       ! Ensure the grid is from -1 to 1, internally.
       x0 = minval(x)
       x1 = maxval(x)
       grid = 2d0*(x-x0)/(x1-x0)-1d0
       
       backward_transform_matrix(:,1) = 0.5d0
       Do j=2,N
          tmp = Dble(j-1)*acos(grid)
          backward_transform_matrix(:,j) = cos(tmp)
       End Do
    End If
  End Subroutine Build_Backward_Transform

  ! Perform the forward transform, using Lapack to compute the solution via dgetrs
  Subroutine Forward_Chebyshev_Transform(x,f,cf)
    Implicit None
    Real*8, Intent(In) :: x(:), f(:,:)
    Real*8, Intent(InOut) :: cf(:,:)
    Integer :: N ,K, info

    If (.not. Allocated(forward_transform_matrix)) Then
       Call Build_Forward_Transform(x)
    End If
    
    cf = f
    N = size(x)
    K = size(f)/N ! Assumes cf and f are the same size
    Call dgetrs('N',N,K,forward_transform_matrix,N,pivots,cf,N,info)
    ! Trap error
    if (info<0) Then
       Print*, 'Illegal value to detrs in forward chevyshev transform at ', info
    End if
  End Subroutine Forward_Chebyshev_Transform

  ! Perform the backward transform, using Lapack via dgemm.
  Subroutine Backward_Chebyshev_Transform(x,cf,f)
    Implicit None
    Real*8, Intent(In) :: x(:), cf(:,:)
    Real*8, Intent(InOut) :: f(:,:)
    Integer :: N ,K

    If (.not. Allocated(backward_transform_matrix)) Then
       Call Build_Forward_Transform(x)
    End If
    
    N = size(x)
    K = size(cf)/N ! Assumes f and cf are the same size
    Call dgemm('N','N',N,K,N,1d0,backward_transform_matrix,N,cf,N,0d0,f,N)
  End Subroutine Backward_Chebyshev_Transform
  
  ! Computes the linking coefficient between aj and uk
  ! In the product between two Ultraspherical polynomial
  ! series of order lambda, e.g. a(x) u(x).
  ! See Olver and Townsend 2013, equation 3.9.
  ! This is constructed to make sure each term is of order 1, to minimize floating point problems.
  ! This could be vectorized over s and k in the next routine.
  Subroutine Multiplication_Coefficient(s,j,k,lambda,csl)
    Implicit None
    Integer, Intent(In) :: s, j, k, lambda
    Real*8 :: csl
    Integer :: t
    Real*8 :: t1, t2, t3, t4, t5
    
    t1 = Dble(j + k + lambda)
    t2 = t1 - Dble(s)
    t1 = t1 - Dble(2*s)
    
    !First term
    csl = t1/t2
    
    t2 = Dble(1)
    Do t=0,s-1
       t2 = t2*Dble(lambda + t)/Dble(1 + t)
    End Do
    
    !Second term
    csl = csl*t2
    
    t3 = Dble(1)
    Do t=0,j-s-1
       t3 = t3*Dble(lambda + t)/Dble(1 + t)
    End Do
    
    !Third term
    csl = csl*t3
    
    t4 = Dble(1)
    Do t=0,s-1
       t1 = Dble(j + k - 2*s + t)
       t2 = t1 + Dble(lambda)
       t1 = t1 + Dble(2*lambda)
       t1 = t1/t2
       t4 = t4*t1
    End Do
    
    !Fourth term
    csl = csl*t4
    
    t5 = Dble(1)
    Do t=0,j-s-1
       t1 = Dble(k - s + t)
       t2 = t1 + Dble(lambda)
       t1 = t1 + Dble(1)
       t1 = t1/t2
       t5 = t5*t1
    End Do
    
    !Fifth term
    csl = csl*t5
  End Subroutine Multiplication_Coefficient

  ! Build the matrix that represents multiplication by a truncated polynomials series
  ! of a scalar field in spectral space. See Olver and Townsend 2013, equation 3.7-3.9.
  Subroutine Build_Multiplication_Matrix(coefs,Ncoefs,N,lambda,matrix)
    Implicit None
    Real*8, Intent(In) :: coefs(:)
    Integer, Intent(In) :: Ncoefs, N, lambda
    Integer :: smin, j, k, s, imin(1), imax(1)
    Integer, Dimension(0:N-1) :: inds, sinds, all_inds
    Real*8 :: mjk, csl
    Real*8, Allocatable :: coefs_copy(:)
    Real*8, Intent(InOut) :: matrix(:,:) 
    
    ! Can save memory here with better choices of data, row, column sizes through heuristics.
    ! For instance, the bandwidth should be about Ncoefs, so maybe 4*Ncoefs*N would be good...
    Allocate(coefs_copy(0:Ncoefs-1))
    coefs_copy(0:Ncoefs-1) = coefs(1:Ncoefs)

    matrix = 0
    
    Do j=0,N-1
       all_inds(j) = j
    End Do
    
    Do j=0,N-1
       Do k=0,N-1
          mjk = 0
          smin = max(0,k-j)
          sinds = -1
          sinds(smin:k) = all_inds(smin:k)
          inds = -1
          inds(smin:k) = 2*all_inds(smin:k)+j-k
          ! Find minimum and maximum compatible indices.
          ! This is to minimize computation and maximize sparsity.
          imin = minloc(inds,mask=(inds.ge.0))
          imax = maxloc(inds,mask=((inds.lt.Ncoefs).and.(inds.ge.0)))
          ! This could be vectorized if necessary.
          If (imin(1).le.imax(1)) Then
             Do s=imin(1),imax(1)
                Call Multiplication_Coefficient(sinds(s),k,inds(s),lambda,csl)
                mjk = mjk + csl*coefs_copy(inds(s))
             End Do
             matrix(j,k) = mjk
          End If
       End Do
    End Do
    
    Deallocate(coefs_copy)
  End Subroutine Build_Multiplication_Matrix

  ! Construct the derivative matrix of order lambda and size N x N, projects T_n -> C^lambda
  Subroutine Build_Derivative_Matrix(N,lambda,matrix)
    Implicit None
    Integer, Intent(In) :: N, lambda
    Real*8, Intent(InOut) :: matrix(:,:)
    Integer :: i
    
    matrix = 0

    Do i=0,N-lambda-1
       matrix(i+1,i+lambda+1) = Dble(i + lambda)
    End Do

    matrix = Dble(2)**(lambda-1)*gamma(Dble(lambda))*matrix
  End Subroutine Build_Derivative_Matrix
  
  ! Construct the matrices that convert coefficients expressed over C^lambda to C^{lambda +1}
  ! S0 converts coefficients of a field a(x) from T_n  to C^1
  ! S(lambda) converts coefficients of a field a(x) from C^lambda to C^{lambda +1}  
  Subroutine Build_Conversion_Matrix(N,lambda,matrix)
    Implicit None
    Integer, Intent(In) :: N, lambda
    Integer :: i, j
    Real*8, Intent(InOut) :: matrix(:,:)

    matrix = 0
    If (lambda.eq.0) Then
       ! Diagonal
       matrix(1,1) = 1 
       Do i=2, N
          matrix(i,i) = 0.5
       End Do
       ! Super-diagonal
       Do i=1, N-2
          matrix(i,i+2) = -0.5
       End Do
    Else
       ! Diagonal
       matrix(1,1) = 1 
       Do i=2, N
          matrix(i,i) = Dble(lambda)/Dble(lambda+i-1)
       End Do
       ! Super-diagonal
       Do i=1, N-2
          matrix(i,i+2) = -Dble(lambda)/Dble(lambda+i+1)
       End Do
    End If
  End Subroutine Build_Conversion_Matrix

  ! This constructs the sparse representation of Dirichlet and Neumann boundary conditions as size N x N.
  ! K is the order of the boundary conditions (e.g. K=4 for a fourth-order ODE).
  ! Here, conditions(1:K/2) specifies the lower boundary, and conditions(K/2+1:K) the upper boundary.
  ! Assuming balanced boundary conditions anyway...
  Subroutine Build_Boundary_Matrix(N,K,conditions,matrix)
    Implicit None
    Integer, Intent(In) :: N, K, conditions(:)
    Integer :: i, j, o
    Real*8, Intent(InOut) :: matrix(:,:)
    Real*8 :: t
    
    matrix = 0
    ! Lower boundary
    Do j=1,K/2
       ! Dirichlet
       If (conditions(j).eq.0) Then
          Do i=0,N-1
             matrix(N-K+j,i+1) = Dble(1 - 2*mod(i,2))
          End Do
       ! Neumann   
       ElseIf (conditions(j).eq.1) Then
          Do i=1,N-1
             matrix(N-K+j,i+1) = Dble((1 - 2*mod(i,2))*i*(i+1))
          End Do
       ! Higher-order p>1   
       Else
          Do i=conditions(j),N-1
             t = 1d0
             Do o=0,conditions(j)-1
                t = t*Dble(i**2-o**2)/Dble(2*o+1)
             End Do
             matrix(N-K+j,i+1) = Dble(1 - 2*(mod(i,2)))*t
          End Do
          matrix(N-K+j,:) = matrix(N-K+j,:)*(-1d0)**conditions(j)
       End If
    End Do

    ! Upper boundary
    Do j=K/2+1,K
       ! Dirichlet
       If (conditions(j).eq.0) Then
          matrix(N-K+j,:) = 1
       ! Neumann   
       ElseIf (conditions(j).eq.1) Then
          Do i=1,N-1
             matrix(N-K+j,i+1) = i*(i+1)
          End Do
       ! Higher-order p>1   
       Else
          Do i=conditions(j),N-1
             t = 1d0
             Do o=0,conditions(j)-1
                t = t*Dble(i**2-o**2)/Dble(2*o+1)
             End Do
             matrix(N-K+j,i+1) = t
          End Do
       End If
       
    End Do
  End Subroutine Build_Boundary_Matrix
  
End Module Ultraspherical_Polynomials_NoSparse
