program main
    !*********************************************************************************************!
    !* Esta rotina encontra o polin�omio que interpola os pontos dados fornecendo os coeficientes *!
    !* num vetor.                                                                                *!
    !*                                                                                           *!
    !* VARIAVEIS:                                                                                *!
    !* ----------                                                                                *!
    !* nnos: N�umero de n�os a serem interpolados                                                  *!
    !* p,q: Vetores arbitr�arios para definir o produto de polin�omios na subrotina ProdPoli       *!
    !* w: Produto dos polin�omios p e q                                                           *!
    !* nos: Vetor contendo os n�os                                                                *!
    !* fVal: Vetor contendo os valores dos n�os                                                   *!
    !* DifDiv: Matriz que cont�em as Diferen�cas divididas dos n�os                                 *!
    !* Lag: Matriz que cont�em os mon�omios do m�etodo de Newton (Diferen�cas divididas)             *!
    !* Interpol: Matriz contendo em cada coluna os termos da interpola��cao de Newton              *!
    !* px: Vetor contendo os coeficientes do polin�omio resultante                                *!
    !*                                                                                           *!
    !* SUBROTINAS:                                                                               *!
    !* -----------                                                                               *!
    !* Dif: Calcula as diferen�cas divididas dos n�os e as armazena na matriz DifDiv. A primeira   *!
    !*      coluna de DifDiv cont�em as zero diferen�cas divididas, a segunda coluna cont�em as pri-*!
    !*      meiras diferen�cas divididas, e assim sucessivamente.                                 *!
    !* ProdPoli: Calcula o produto de dois polin�omios cujos coeficientes s�ao armazenados nos     *!
    !*           vetores p e q. A primeira componente destes vetores �e o termo constante, a se-  *!
    !*           gunda o coeficiente de x, a terceira o coeficiente de x^3, etc. Os coeficientes *!
    !*           do produto p*q s�ao armazenados no vetor w.                                      *!
    !*********************************************************************************************!
    implicit none
    integer :: nnos = 5,j
    real :: soma
    real, allocatable :: p(:), q(:), w(:), nos(:), fVal(:), px(:), soma1(:)
    real, allocatable :: DifDiv(:,:), Lag(:,:), Interpol(:,:)
    allocate(nos(nnos),fVal(nnos),w(nnos),px(nnos),soma1(nnos))
    allocate(DifDiv(nnos,nnos),Lag(nnos,nnos),Interpol(nnos,nnos))


    DifDiv = 0
    Lag = 0
    soma1 = 0
    nos = (/ 1., 2., 4., 6., 8. /)
    fVal = (/ 1., 3., 6., 10., 12. /)
    DifDiv(:,1) = fVal
    Lag(1,1) = 1

    do j = 2,nnos
        Lag(1,j) = -nos(j-1)
        Lag(2,j) = 1
        call ProdPoli(Lag(:,j-1),Lag(:,j),w)
        Lag(:,j) = w
    end do

    call Dif(nos,fVal,DifDiv)

    do j = 1,nnos
        Interpol(:,j) = DifDiv(1,j)*Lag(:,j)
    end do

    do j = 1,nnos
        soma1 = soma1 + Interpol(:,j)
    end do

    px = soma1
    write(*,*)px

    contains

    subroutine Dif(nos,fVal,DifDiv)
        real, allocatable :: nos(:), fVal(:)
        real, allocatable :: DifDiv(:,:)
        integer :: i,j
        do j = 2,nnos
            do i = 1,nnos-(j-1)
                DifDiv(i,j) = (DifDiv(i+1,j-1)-DifDiv(i,j-1))/(nos(i+j-1) - nos(i))
            end do
        end do
    end subroutine

    subroutine ProdPoli(p,q,w)
        real :: soma
        real :: p(:), q(:), w(:)
        integer :: k, jmax, jmin, j

        do k = 1,size(p)+size(q)-1
            soma = 0
            jmax = max(1,k+1-size(p))
            jmin = min(k,size(q))
            do j = jmax,jmin
                soma = soma + p(j)*q(k-j+1)
            end do
            w(k) = soma
        end do
    end subroutine ProdPoli

end program main
