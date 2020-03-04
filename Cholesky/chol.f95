program cholesky

    implicit none
    integer, parameter :: n = 4 !LINHAS
    real*8, dimension(n,n) :: L !ORDEN DA MATRIZ TRIANGULAR L
    real*8, dimension(n,n) :: A !ORDEN DA MATRIZ INICIAL
    real*8 :: soma1, soma2, soma3
    integer :: i, j, k

    !MATRIZ INICIAL
    A(1,:) = (/ 6., 2., 1., -1. /)
    A(2,:) = (/ 2., 4., 1., 0. /)
    A(3,:) = (/ 1., 1., 4., -1. /)
    A(4,:) = (/ -1., 0., -1., 3. /)

    L(1,1) = sqrt(A(1,1))

    do i = 2,n
        L(i,1) = A(i,1)/L(1,1)
    end do

    do k = 2,n
        soma1 = 0
        do j = 1,k-1
            soma1 = soma1 + L(k,j)**2
        end do
        L(k,k) = sqrt(A(k,k)-soma1)
        do i = k+1,n
            soma2 = 0
            do j = 1,k-1
                soma2 = soma2 + L(k,j)*L(i,j)
            end do
            L(i,k) = (A(i,k)-soma2)/L(k,k)
        end do
    end do

    print*, "L = "
    do i = 1,n
        write(*,*)L(i,:)
    end do
    
end program cholesky
