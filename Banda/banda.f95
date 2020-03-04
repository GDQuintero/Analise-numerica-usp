program banda

    implicit none
    integer, parameter :: n = 6 !LINHAS E COLUNAS (MATRIZ QUADRADA)
    integer, parameter :: bn = 5 !TAMANHO DA BANDA OU NUMERO DE DIAGONAIS
    real*8, dimension(n,n) :: A
    real*8, dimension(n,bn) :: B = 0.d0
    real*8 :: m
    integer :: i, j, k

    !MATRIZ DE BANDA
    A(1,:) = (/ 5.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0 /)
    A(2,:) = (/ 1.d0, 5.d0, 1.d0, 1.d0, 0.d0, 0.d0 /)
    A(3,:) = (/ 1.d0, 1.d0, 5.d0, 1.d0, 1.d0, 0.d0 /)
    A(4,:) = (/ 0.d0, 1.d0, 1.d0, 5.d0, 1.d0, 1.d0 /)
    A(5,:) = (/ 0.d0, 0.d0, 1.d0, 1.d0, 5.d0, 1.d0 /)
    A(6,:) = (/ 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, 5.d0 /)


    !MATRIZ AUXILIAR

    !DIAGONAL PRINCIPAL DA MATRIZ A
    do i = 1,n
        B(i,(bn+1)/2) = A(i,i)
    end do

    !CODIAGONAIS DA MATRIZ A
    do k = 1,(bn-1)/2
        do i = 1,(n-k)
            B(i+k,-k+(bn+1)/2) = A(i+k,i)
            B(i,k+(bn+1)/2) = A(i,i+k)
        end do
    end do

    !IMPRIMINDO MATRIZ AUXILIAR
    print*
    print*, "Matrix auxiliar inicial: "
    print*
    do i = 1,n
        write(*,*)B(i,:)
    end do

    !APLICANDO ELIMINACAO GAUSSIANA NA MATRIZ AUXILIAR
    do k = 1,n-1
        do i = k+1,k+(bn-1)/2
            m = B(i,k-i+(bn+1)/2)/B(k,(bn+1)/2)
            do j = k+1,i-1+(bn-1)/2
                B(i,j-i+(bn+1)/2) = B(i,j-i+(bn+1)/2) - m*B(k,j-k+(bn+1)/2)
            end do
            B(i,k-i+(bn+1)/2) = 0
        end do
    end do

    print*
    print*,"Matriz auxiliar apos aplicar eliminacao Gaussiana: "
    print*

    do i = 1,n
        write(*,*)B(i,:)
    end do


end program banda
