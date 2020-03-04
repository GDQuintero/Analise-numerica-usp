PROGRAM LU
    integer, parameter :: n = 4
    real*8, dimension(n,n) :: A, L=0, U=0
    real*8 :: soma1, soma2

    A(1,:) = (/ 1, 1, 0, 3 /)
    A(2,:) = (/ 2, 1, -1, 1 /)
    A(3,:) = (/ 3, -1, -1, 2 /)
    A(4,:) = (/ -1, 2, 3, -1 /)

    U(1,:) = A(1,:)

    do i = 1,n
        L(i,i) = 1
    end do

    do i = 2,n
        L(i,1) = A(i,1)/U(1,1)
    end do

    do k = 2,n
        do j = k,n
            soma1 = 0
            do s = 1,k-1
                soma1 = soma1 + L(k,s)*U(s,j)
            end do
            U(k,j) = A(k,j) - soma1
        end do
        do i = k+1,n
            soma2 = 0
            do s = 1,k-1
                soma2 = soma2 + L(i,s)*U(s,k)
            end do
            L(i,k) = (A(i,k) - soma2)/U(k,k)
        end do
    end do

    print*, "L = "
    do i = 1,n
        write(*,*)L(i,:)
    end do

    print*

    print*, "U = "
    do i = 1,n
        write(*,*)U(i,:)
    end do
END PROGRAM LU
