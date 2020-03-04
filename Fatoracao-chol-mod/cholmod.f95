PROGRAM cholmod

    implicit none
    integer, parameter :: n = 3
    real*8, dimension(n,n) :: A
    real*8, dimension(n,n) :: D=0, L=0
    real*8 :: soma1, soma2
    integer :: i, j, k

    A(1,:) = (/ 4., 1., 1. /)
    A(2,:) = (/ 1., 2., -1. /)
    A(3,:) = (/ 1., -1., 3. /)

    D(1,1) = A(1,1)

    do i = 1,n
        L(i,i) = 1
    end do

    do i = 2,n
        L(i,1) = A(i,1)/D(1,1)
    end do

    do k = 2,n
        soma1 = 0
        do i = 1,k-1
            soma1 = soma1 + (L(k,i)**2)*D(i,i)
        end do
        D(k,k) = A(k,k) - soma1
        do j = k+1,n
            soma2 = 0
            do i = 1,k-1
                soma2 = soma2 + L(k,i)*L(j,i)*D(i,i)
            end do
            L(j,k) = (A(j,k) - soma2)/D(k,k)
        end do
    end do

    print*, "L = "
    do i = 1,n
        write(*,*)L(i,:)
    end do
    print*
    print*, "D = "
    do i = 1,n
        write(*,*)D(i,:)
    end do
END PROGRAM cholmod
