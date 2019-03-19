Program Driver_HH

  ! load module and specific which functions and variables to use in the driver
  use LinAl, only: mat, msize, nsize, readMat, printMat,vectorTwoNorm,findTrace, GEPivot, &
        backsub, LUDecomp, LUBacksub, HHbuildA, HHQR

  implicit none
  
  character(len=100) :: myFileName
  real(kind =8), dimension(:,:), allocatable:: Amat, Bmat, Qmat, Emat1, Emat2, EmatQ, &
          Imat, xmat, ymat, As, Bs
  real(kind = 8) :: p_max, k, trace
  integer :: i,j, Ndegree
  
  myFileName = 'atkinson.dat'

  ! degree of polynomial fit
  Ndegree = 5

  ! load A-matrix
  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  allocate(Amat(msize, Ndegree+1))
  allocate(As(msize, Ndegree+1))
  allocate(Imat(msize,msize))
  allocate(Qmat(msize, msize))
  allocate(EmatQ(msize, msize))
  allocate(Bmat(msize, 1))
  allocate(Bs(msize,1))
  allocate(Emat1(msize, Ndegree+1))
  allocate(Emat2(msize, 1))
  allocate(xmat(Ndegree+1, 1))
  allocate(ymat(msize, 1))

  ! Always initialize with zeros
  mat = 0.0
  call readMat(myFileName) ! call the function and load the matrix
  call HHbuildA(mat, msize, Ndegree, Amat) !constructed A-matiix
  As = Amat
  Bmat = reshape(mat(:, 2), (/msize,1/)) !constructed B-matrix (2nd column of the data file)
  Bs = Bmat

  ! print the matrix in readable format
  write(*,*) "Input A-matrix"
  call printMat(Amat)
  write(*,*) "--------------------------------------------------------------"
  write(*,*) "Input B-matrix"
  call printMat(Bmat)
  write(*,*) "--------------------------------------------------------------"

  ! implement HH QR decomposion
  call HHQR(Amat, msize, Ndegree+1, Qmat)
  write(*,*) "Reduced R-matrix"
  call printMat(Amat)
  write(*,*) "--------------------------------------------------------------"

  ! calculate the error of QR decomposition
  Emat1 = matmul(Qmat,Amat)
  Emat1 = As - Emat1
  write(*,*) "Output Error-matrix : A - QR"
  call printMat(Emat1)
  write(*,*) "--------------------------------------------------------------"

  write(*,*) "Frobenius norm of A-QR"
  write(*,*) norm2(Emat1)
  write(*,*) "--------------------------------------------------------------"

  ! calculate the error of Q, Q^TQ - I
  Imat = 0.0
  do i = 1, msize
    Imat(i,i) =1
  end do
  EmatQ = matmul(transpose(Qmat), Qmat)
  EmatQ = EmatQ - Imat

  write(*,*) "Output Error-matrix : Q^TQ - I "
  call printMat(EmatQ)
  write(*,*) "--------------------------------------------------------------"

  write(*,*) "Frobenius norm of Q^TQ-I "
  write(*,*) norm2(EmatQ)
  write(*,*) "--------------------------------------------------------------"

  ! solve the least square problem Ax = b >> QRx = b >> Rx = Q^Tb >> Rx = y ::
  ymat = matmul(transpose(Qmat), bmat) ! y = Q^Tb
  xmat = 0.0
  write(*,*) "Least square problem y-matrix "
  call printMat(ymat)
  write(*,*) "--------------------------------------------------------------"

  ! solve Rx = y:
  call backsub(Amat, ymat, Ndegree+1, 1, xmat)
  write(*,*) "LSP solution x-matrix "
  call printMat(xmat)
  write(*,*) "--------------------------------------------------------------"

  ! verify x-mat by calculating the error: bmat - Ax
  Emat2 = matmul(As, xmat)
  Emat2 = bmat - Emat2

  write(*,*) "Error matrix for b - Ax "
  call printMat(Emat2)
  write(*,*) "--------------------------------------------------------------"

  write(*,*) " 2norm of b - Ax"
  write(*,*) norm2(Emat2)
  write(*,*) "--------------------------------------------------------------"

  deallocate(mat)
  deallocate(Amat)
  deallocate(Bmat)
  deallocate(Qmat)
  deallocate(EmatQ)
  deallocate(As)
  deallocate(Bs)
  deallocate(Emat1)
  deallocate(xmat)
  deallocate(ymat)
  deallocate(Emat2)



End Program Driver_HH
