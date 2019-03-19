module LinAl
  implicit none
  
  integer, save :: msize, nsize
  real, dimension(:,:), allocatable, save :: mat

contains

  !********************************************************
  subroutine printMat(mat)
    ! print out matrix in a human readable format

    implicit none
    real, dimension(:,:):: mat
    integer, dimension(2) :: mn
    integer :: i,j

    mn = shape(mat)
    do i = 1, mn(1)
      write(*,*) (mat(i,j) , j = 1, mn(2) )
    end do

  end subroutine printMat

  !********************************************************

  subroutine readMat(filename)

    implicit none
    character(len=*) :: filename
    integer :: i,j, k

    ! Reads a file containing the matrix A 
    ! Sample file:
    ! 3  4
    ! 1.2     1.3     -1.4    3.31
    ! 31.1    0.1     5.411   -1.23
    ! -5.4    7.42    10      -17.4
    ! Note that the first 2 lines are the matrix dimensions, 
    ! then the next msize lines are the matrix entries
    ! Note that entries must be separated by a tab.


    open(10,file=filename)

    ! Read the matrix dimensions -- here, this part is a dummy part and don't do anything
    read(10,*) i,j

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)

  end subroutine readMat

  !********************************************************

  subroutine findTrace(mat,msize,tr)
    ! Find the trace of the a matrix: sum of the diagonal
    ! input: any squared matrix; msize= # of row for matrix
    ! output: trace of the matrix
    implicit none

    real(kind = 8), dimension(:,:), intent(in):: mat
    integer,  intent(in):: msize
    integer :: i
    real(kind =8), intent(out) :: tr

    tr = 0.0
    do i =1, msize
      tr = tr + mat(i,i)
    enddo

  end subroutine findTrace

  !********************************************************

  subroutine vectorTwoNorm(Emat, m, n, norm2)
    ! calculate 2norm of each column of the input error matrix
    ! input: Emat = error matrix ; m = # of row of matrix; n = # of column of matrix
    ! output: norm2 : vector contains the two norm of each column

    implicit none

    integer, intent(in)::m,n
    integer :: j
    real(kind = 8), dimension(n):: vsum
    real(kind = 8), dimension(m,n), intent(in) :: Emat
    real(kind = 8), dimension(n), intent(out):: norm2

    vsum = 0.0
    do j =1, n
      vsum(:) = vsum(:) + Emat(:,j)**2
    end do
    norm2(:) = sqrt(vsum(:))

  endsubroutine vectorTwoNorm

  !********************************************************

  subroutine  GEPivot(Amat, Bmat, msize)
    ! implementing Gaussian elimination with partial pivoting for linear eqn Ax=B
    ! input: Amat: matrix A; Bmat: matrix B; msize: row size of A & B
    ! output:   Upper triangular matrix U that stored in A, Updated B matrixed.

    implicit none

    real(kind =8), dimension(:,:), intent(inout) :: Amat, Bmat
    real(kind=8) :: a_ratio
    integer, dimension(2) :: mn
    integer :: i, j, k, msize
    real(kind=8), dimension(msize) :: temp

    mn = shape(Amat)
    do j = 1, mn(1) -1
      ! implement partial pivoting
      k = maxloc(Amat(:, j), dim=1) !maxloc search the maximium value of a array and return the index
      if (k > j) then
        ! interchange rows for A-matrix
        temp = Amat(k,:) !storage k-row vector in a dummy variable
        Amat(k,:) = Amat(j,:) !replace k-row by j-row
        Amat(j,:) = temp ! replace j-row by the storaged k-row
        ! interchange rows for B-matrix
        temp = Bmat(k,:) !storage k-row vector in a dummy variable
        Bmat(k,:) = Bmat(j,:) !replace k-row by j-row
        Bmat(j,:) = temp ! replace j-row by the storaged k-row
      end if
      ! check if A-matrix is singular or not
      if (Amat(j,j) == 0) then
        stop ! A-matrix is singular
      endif
      ! implement Gaussian Elimination
      do i = j+1, mn(1)
        ! eliminate the lower triangular element of A-matrix into zero
        a_ratio = Amat(i,j) / Amat(j,j)
        Amat(i, :) = Amat(i, :) - Amat(j, :) * a_ratio ! A(i,:) support vector operation
        Bmat(i, :) = Bmat(i, :) - Bmat(j, :) * a_ratio
      enddo
    enddo

  end subroutine GEPivot

  !********************************************************

  subroutine backsub (Umat, Bmat, Bmsize, Bnsize, Xmat)
    ! perform GE back-substitution to solve Ax=B
    ! input: Umat: triangulized A-matrix, Bmat: B-matrix after GE;
    !             Bmsize, # of row of B-matrix; Bnsize: # of column of B-matrix
    ! output Xmat: Solution matrix

    implicit none

    real(kind=8), dimension(:,:), intent(in):: Umat, Bmat
    real(kind=8), dimension(Bmsize, Bnsize),intent(out) :: Xmat
    integer, intent(in) :: Bmsize, Bnsize
    integer:: i, j
    real(kind=8), dimension(Bnsize):: sum

    ! first column
    Xmat(Bmsize, :) = Bmat(Bmsize, :)/Umat(Bmsize, Bmsize) !calculate the solution of last row

    do i = Bmsize -1, 1, -1 !row iteration start from the bottom
      sum = 0.0
      do j = i+1, Bmsize ! column iteration
        sum(:) = sum(:) + Umat(i,j) * Xmat(j,:) !sum up all the values from the solved terms
      end do
      Xmat(i, :) = (Bmat(i,:) - sum(:))/Umat(i,i)
    end do



  end subroutine backsub

  !********************************************************

  subroutine LUDecomp (Amat, msize, S)
    ! perform LU decomposition on input matrix A
    ! input: Amat: input A-matrix; msize: first dimension of A-matrix (# of row)
    ! output: S: permutation vector records the pvioting

    implicit none

    real(kind =8), dimension(:,:), intent(inout):: Amat
    real(kind=8), dimension(msize):: S, temp
    real(kind=8) :: dummy
    integer, intent(in) :: msize
    integer:: i,j,k


    do j = 1, msize
      S(j) = j !create permutation vector, S
    enddo

    do j =1, msize
      ! implement partial pivoting
      k = maxloc(Amat(:, j), dim=1) !maxloc search the maximium value of a array and return the index
      if (k > j) then
        ! interchange rows for A-matrix
        temp = Amat(k,:) !storage k-row vector in a dummy variable
        Amat(k,:) = Amat(j,:) !replace k-row by j-row
        Amat(j,:) = temp ! replace j-row by the storaged k-row
        ! interchange rows for S-vector
        dummy = S(k) !storage k element in a dummy variable
        S(k) = S(j) !replace k-element by j-element
        S(j) = dummy ! replace j-element by the storaged k-element
      end if
      if (Amat(j,j) == 0) then
        stop ! stop if A is singular
      end if
      do i = j+1, msize
        Amat(i,j) = Amat(i,j)/Amat(j,j) ! update A and save l(i,j) into A(i,j)
        do k = j+1, msize
          Amat(i,k) = Amat(i,k) - Amat(i,j) * Amat(j,k) !update A by results of elimination (U)
        end do
      end do
    end do

    ! compact way of LU decomp without pivoting
    ! source: https://rosettacode.org/wiki/LU_decomposition#Fortran
!    do k = 1,msize-1
!      Amat(k+1,k) = Amat(k+1,k) / Amat(k,k)
!      forall (j=k+1:msize) Amat(k+1,j) = Amat(k+1,j) - Amat(k+1,k) * Amat(k,j)
!    end do

  end subroutine LUDecomp

  !********************************************************


  subroutine LUBacksub (Amat, Bmat, m, n, S, Ymat, Xmat)
    ! perform LU forward substitution & backward substitution
    ! input: Amat: decomposed A-matrix that contain U & L matrices; Bmat: original B-matrix;
    !               m: # of row of B-matrix; n: # of row of B-matrix; S: permutation vector;
    ! output: Ymat: Solution of forward sub. in Ly = b;
    !         Xmat: solution of backward sub. in Ux = y.

    implicit none

    real(kind = 8), dimension(:,:), intent(in) :: Amat, Bmat
    real(kind=8), dimension(:), intent(in) :: S
    real(kind=8), dimension(:,:), intent(out) :: Ymat, Xmat
    real(kind=8), dimension(n) :: sum
    integer, intent(in) :: m, n
    integer :: i,j,k,a

    ! create permutated B-matrix and save in Y-vector
    do j=1,m
      a = S(j) ! convert real value in S into integer for indicies
      Ymat(j,:) = Bmat(a, :)
    end do

    ! Forward Subsitution: solve y in Ly = Pb
    do j =1, m -1
      do i = j+1, m
        Ymat(i, :) = Ymat(i, :) - Ymat(j,:) * Amat(i,j)
      end do
    end do

    ! Backward Subsitution: solve x in Ux = y
    do i = m,1,-1
      sum =0.0
      do k= i+1, m
        sum(:) = sum(:) + Amat(i,k)*Xmat(k, :)
      end do
      Xmat(i,:) = (Ymat(i,:) - sum(:))/Amat(i,i)
    end do

  end subroutine LUBacksub

  !********************************************************

  subroutine HHbuildA(data, m, n, Amat)
    !countruct a polynimal A-matrix based on the input data file for
    !     Householder QR decomp. fit
    !   input: data = given data file; m = rowsize of A-matrix; n = degree of polynomial;
    !   output: Amat = counstructed A-matrix

    implicit none

    real(kind =8), dimension(:,:), intent(in) :: data
    real(kind =8), dimension(m, n+1), intent(out) :: Amat
    integer, intent(in):: m,n
    integer:: i,j

    Amat = 0.0
    Amat(:,1) = 1
    do j = 2, n+1
      Amat(:,j) = data(:,1)**(j-1) !replace each column of A with polynomial
    end do
  end subroutine HHbuildA

  !********************************************************

  subroutine HHQR (Amat, m, n, Qmat)
    ! perform Householder QR decomposition of the input matrix.
    ! input: Amat= input matrix (mxn); m = row size of A-matrix
    !        n = colunm size of A-matrix;
    ! output: Amat = reduced A-matrix which is R-matrix; Qmat: Q-matrix with HH vectors

    implicit none

    real(kind =8), dimension(:,:), intent(inout) :: Amat
    real(kind =8), dimension(m,m), intent(out):: Qmat
    real(kind =8), dimension(m,m) ::Imat, Hmat
    real(kind =8), dimension(m):: v,a
    real(kind =8):: s
    real(kind =8)::normA, normv
    integer, intent(in) :: m, n
    integer::i , j

     ! create an identity matrix of size mxm
    Imat = 0.0
    do i = 1,m
      Imat(i,i) = 1
    end do
    Qmat = Imat

    do j = 1,n ! column iterator

      a = Amat(:,j)

      ! check j is already more than 1 step, if yes, insert 0!!!
      if (j-1 .NE. 0) then
        ! insert 0 into a-vector to make sure the previou row of A-matrix stay unchaged
        do i = 1,j-1
          a(i) = 0
        end do
      end if

      normA = 0.0
      ! calculate the norm of j-column of A-matrix
      do i = 1, m
        normA = normA + a(i) ** 2
      end do
      normA = sqrt(normA)

      ! determine the sign
      s = sign(normA, Amat(j,j))

      ! construct v-vector
      v = a
      v(j) = Amat(j,j) + s

      ! calculate the norm of v
      normv = 0.0
      do i = 1, m
        normv = normv + v(i) ** 2
      end do
      normv = sqrt(normv)
      v = v / normv

      ! update A-matrix for QR-reduction
      ! notice: "reshape" function change the dimension of the vector/matrix
      Amat(:,j:n) = Amat(:,j:n) - 2*matmul(reshape(v, (/m, 1/)), matmul(reshape(v,(/1,m/)),Amat(:,j:n)))

!      write(*,*) "-------------"
!      call printMat(Amat)

      ! construct Q-matrix
      ! save H-matrix for each step
      Hmat = Imat -  2* matmul(reshape(v, (/m, 1/)), reshape(v,(/1,m/)))
      ! accumulate H-matrix from each step
      Qmat = matmul(Hmat, Qmat)
    end do

    Qmat = transpose(Qmat)

  end subroutine HHQR

  !********************************************************

  subroutine HH_Tridiag (Amat, m, n)
    ! only works for square matrix since the aim is to calculate eigenvenvalue.
    ! perform Similar transformation on input matrix by HH project;
    ! note: for symmetry matrix, it will produce tridiagonal matrix;  otherwise it will gives rise Hessenberg in general
    !
    ! input: Amat = squared matrix; m = n: row and column size of A matrix;
    ! output: Amat = transformed A-matrix that is either tridiagonal or Hessenberg.

    implicit none

    real(kind = 8), dimension(:,:), intent(inout) :: Amat
    real(kind = 8), dimension(m) :: a, v
    real(kind = 8) :: s, normA, normv
    integer :: m, n , i, j , k

    do k = 1, m-2

      a = 0.0 ! always intiate with 0 value
      v = 0.0
      j = k + 1! dummy variable

      a(1: m-k) = Amat(j:m, 1) ! pick the 1st column of A-matrix, since we are subseting A each time

      ! calculate the 2norm of vector a
      normA = 0.0
      do i = 1, m-k
        normA = normA + a(i) ** 2
      end do
      normA = sqrt(normA)

      ! construct vector v
      v = a
      v(1) = a(1) + sign(normA, a(1))
      ! since we perfrom subseting, so we always modify the first element of v
      normv = 0.0  ! calculate 2norm of v
      do i = 1, m-k
        normv = normv + v(i) ** 2
      end do
      normv = sqrt(normv)
      v = v / normv

      ! update A-matrix twice through similarity transformation of HH projector

      Amat(j:m,k:m) = Amat(j:m,k:m) - 2*matmul(reshape(v(1:m-k), (/m-k, 1/)),matmul(reshape(v(1:m-k),(/1,m-k/)),Amat(j:m,k:m)))
      Amat(:,j:m) = Amat(:,j:m) - 2*matmul(matmul(Amat(:,j:m),reshape(v(1:m-k), (/m-k,1/))),reshape(v(1:m-k),(/1,m-k/)))

    end do

  end subroutine HH_Tridiag

  !********************************************************

  subroutine QR_Alg (Amat, m, n, shift)
    ! perform QR algorithms for both with and without shift
    ! that calculate the eigenvalues of an input matrix
    ! inputs: Amat= input squared matrix, m=n = row and colunm size of A;
    !         shift (logical) = True > QR algorithm with shift, False > without shift
    ! ouput:  Aout = diagonal matrix with eigenvalues at the digonal elements

    implicit none

    real(kind = 8), dimension(:,:), intent(inout) :: Amat
    real(kind = 8), dimension(m,n) :: Imat, Qmat, Bmat, Aout
    real(kind = 8) :: u, r
    integer, intent(in) :: m, n
    integer :: i, j , k
    logical, intent(in):: shift

    Qmat = 0.0
    Aout = 0.0

    k= 0
    u = 0
    r = 10.e16 ! initiate with a large residue error
    k = 0

    ! construct an identity matrix
    Imat = 0.0
    do i = 1, m
      Imat(i,i) = 1
    end do



    if (shift .eqv. .TRUE.) then
      ! implement QR algorithm with shift + deflation(not working?1)
      do j = m, 2, -1

        ! perform deflation
        Bmat(1:j, 1:j) = Amat(1:j, 1:j)

        do while (k < 10)
          u = Bmat(j,j) ! always pick the last diagonal element as the shift
          Bmat(1:j, 1:j) = Bmat(1:j,1:j) - u*Imat(1:j, 1:j)

          call HHQR(Bmat(1:j, 1:j), j, j, Qmat(1:j, 1:j) )
          Bmat(1:j,1:j) = matmul(Bmat(1:j, 1:j), Qmat(1:j,1:j)) + u*Imat(1:j,1:j)
          r = norm2(abs(Aout - Amat))
          Amat(1:j, 1:j) = Bmat(1:j,1:j)
          Aout(1:j, 1:j) = Bmat(1:j,1:j)
          k = k + 1

        end do


      end do


    else

      ! perform QR algorithm without shift
      do while (r > 1.e-16)
        call HHQR(Amat, m, n, Qmat)
        Amat = matmul(Amat, Qmat)
        r = norm2(abs(Aout - Amat)) ! calculate the residuse
        Aout = Amat

      end do
    end if

  end subroutine QR_Alg

  !********************************************************


  subroutine II_EigenV (Amat, m, n, u, v)
    ! perform inverse iterative algorithm to determine the eigenvector of the known
    !   eigenvalues
    ! inputs: Amat = input suqared matrix; m = n = row/column size of A;
    !         u = known or approximated eigenvalue(close to the actual one)
    ! output: v = eigenvector

    implicit none

    real(kind = 8), dimension(:,:), intent(in):: Amat
    real(kind = 8), dimension(m,m) :: Imat, Bmat
    real(kind = 8), dimension(m,1), intent(out) :: v
    real(kind = 8), dimension(m) :: w
    real(kind = 8), intent(in) :: u
    real(kind = 8) :: r
    integer, intent(in) :: m,n
    integer :: i, j , k

    ! intial guess of eigenvector v and a close eigenvalue guess u
    v(:,1)= (/1, 0, 0, 0 /)

    ! the extra dimension of v is to match the dimension of previous GEPivet algorithm
    k = 0
    r = 1.e16

    ! create the identity matrix in mxm
    Imat = 0.0
    do i = 1, m
      Imat(i, i) = 1.0
    end do

    ! create the shifted A-matrix (A - uI)
    Bmat = Amat - u*Imat

    !do while (r > 1.e-16)
    do while (k < 20)
      v(:,1) = v(:,1) / norm2(reshape(v, (/m, 1/))) !normalized v

      ! constructed (A-uI)^-1 by solving (A-uI)w = v using Gaussian Elimation & backsub
      call GEPivot(Bmat, v, m)
      call backsub(Bmat, v, m, 1, w)
      ! w is the solution and w = (A-uI)^-1 * v
      w = w / norm2(reshape(w, (/m,1/)))

      ! calculate the residuse for every iteration
      r = norm2(abs(w) - abs(v(:, 1)))

      v(:,1) = -w !the sign is flip during the iteration somehow :/
      k = k + 1
    end do


  end subroutine II_EigenV







end module LinAl
