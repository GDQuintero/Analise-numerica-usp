program Projeto
implicit none
    
    integer :: PtosInt = 7, Ordem = 4, NumBanda, PtosMalha = 10000, i, j, l, n
    real(kind=8), parameter :: pi = 4.d0*Atan(1.d0)
    real(kind=8) :: a = 0.d0, b = 1.d0 !CONDICOES DE FRONTEIRA
    real(kind=8) :: x, h, soma, h1
    real(kind=8), allocatable :: Part(:), PartExt(:), QSplines(:,:,:,:), BaseSplines(:,:,:), BCondCont(:,:,:), SistemaLinear(:,:)
    real(kind=8), allocatable :: VetorSist(:), VetorSol(:), Solucao(:), Malha(:), SolExata(:)
    
    allocate(Part(PtosInt+2),PartExt(PtosInt+2*Ordem),BaseSplines(Ordem+PtosInt,Ordem,Ordem),BCondCont(Ordem+PtosInt-2,Ordem,Ordem))
    allocate(SistemaLinear(Ordem+PtosInt-2,Ordem+PtosInt-2),VetorSist(Ordem+PtosInt-2),VetorSol(Ordem+PtosInt-2))
    allocate(Solucao(PtosMalha+1),Malha(PtosMalha+1),SolExata(PtosMalha+1))
    
    h = real(Ordem)/real((PtosInt+1))
    NumBanda = 2*Ordem-1
    Open(Unit = 10, File = "Resultados.txt")
    SistemaLinear = 0.d0; VetorSist = 0.d0; soma = 0.d0
    QSplines = BSplines(Ordem,PtosInt)
    
    
    !================================================================================================
    ! CONDICOES DE CONTORNO
    !================================================================================================
    !BSPLINES CUBICOS
    if (Ordem .eq. 4) then
        BaseSplines= h*QSplines(2,:,:,:)
        BCondCont = BaseSplines(2:Ordem+PtosInt-1,:,:)

        BCondCont(1,1,:) = 0.d0
        BCondCont(1,2,:) = 0.d0
        BCondCont(1,3,:) = SomPoli(BaseSplines(2,3,:),-4.d0*BaseSplines(1,4,:))
        BCondCont(1,4,:) = BaseSplines(2,4,:)

        BCondCont(2,1,:) = 0.d0
        BCondCont(2,2,:) = SomPoli(BaseSplines(3,2,:),-1.d0*BaseSplines(1,4,:))
        BCondCont(2,3,:) = BaseSplines(3,3,:)
        BCondCont(2,4,:) = BaseSplines(3,4,:)

        BCondCont(Ordem+PtosInt-3,3,:) = SomPoli(BaseSplines(Ordem+PtosInt-2,3,:),-1.d0*BaseSplines(Ordem+PtosInt,1,:))
        BCondCont(Ordem+PtosInt-3,4,:) = 0.d0

        BCondCont(Ordem+PtosInt-2,2,:) = SomPoli(BaseSplines(Ordem+PtosInt-1,2,:),-4.d0*BaseSplines(Ordem+PtosInt,1,:))
        BCondCont(Ordem+PtosInt-2,3,:) = 0.d0
        BCondCont(Ordem+PtosInt-2,4,:) = 0.d0
        
    elseif (Ordem .eq. 2) then
        BaseSplines= h*QSplines(1,:,:,:)
        BCondCont = BaseSplines(2:Ordem+PtosInt-1,:,:)
    end if
    
    !================================================================================================
    ! MONTAGEM E SOLUCAO DO SISTEMA LINEAR
    !================================================================================================
    !DIAGONAL DA MATRIZ
    do i = 1, Ordem+PtosInt-2
        SistemaLinear(i,i) = ProdIntL(BCondCont(i,:,:),BCondCont(i,:,:),i,i,Ordem,PtosInt)
    end do
    
    !ELEMENTOS NAO DIAGONAIS
    do l = 1, Ordem-1
        i = 0
        do j = l+1, Ordem+PtosInt-2
            i = i+1
            SistemaLinear(i,j) = ProdIntL(BCondCont(i,:,:),BCondCont(j,:,:),i,j,Ordem,PtosInt)
            SistemaLinear(j,i) = SistemaLinear(i,j)
        end do
    end do
    
    !VETOR DO SISTEMA
    do i = 1, Ordem+PtosInt-2
        VetorSist(i) = ProdIntUsual(BCondCont(i,:,:),i)
    end do

    !SOLUCAO DO SISTEMA

    VetorSol = GaussBanda(SistemaLinear,VetorSist,NumBanda)
    
    !================================================================================================
    ! EXPORTAMOS A SOLUCAO APROXIMADA E A SOLUCAO EXATA NUM ARQUIVO .TXT
    !================================================================================================  
    h1 = 1.d0/(PtosMalha); Malha(1) = 0.d0
    Solucao(1) = u(Malha(1),BCondCont,VetorSol)
    SolExata(1) = Exata(Malha(1))
    if ((a .ne. 0.d0).or.(b .ne. 1.d0)) then
            Solucao(1) = Solucao(1) + a + (b-a)*Malha(1)
        end if
    do i = 2, PtosMalha+1
        Malha(i) = Malha(i-1)+h1
        Solucao(i) = u(Malha(i),BCondCont,VetorSol)
        if ((a .ne. 0.d0).or.(b .ne. 1.d0)) then
            Solucao(i) = Solucao(i) + a + (b-a)*Malha(i)
        end if
        SolExata(i) = Exata(Malha(i))
    end do

    do i = 1, PtosMalha+1
        write(10,*)Malha(i), Solucao(i), SolExata(i)
    end do
    Close(10)
    print*, Norma(Solucao,SolExata)
    
   contains
    
    !================================================================================================
    ! AVALIACAO DA FUNCAO APROXIMADA U
    !================================================================================================
    function u(x,phi,b)
        implicit none
        
        integer :: i,j,ind
        real(kind=8) :: phi(:,:,:), b(:), x, u, h, soma
        real(kind=8), allocatable :: aux(:)
        
        call Particao(Ordem,PtosInt,Part,PartExt)
        allocate(aux(Ordem+PtosInt-2))
        aux = PartExt(2:2*Ordem+PtosInt-1)
       
        if (Ordem .eq. 4) then
            soma = 0.d0
            do i = 1, Ordem+PtosInt-2
                if (aux(i) .le. x .and. x .le. aux(i+Ordem)) then
                    if (aux(i) .le. x .and. x .le. aux(i+1)) then
                        soma = soma + VetorSol(i)*Horner(phi(i,1,:),x,0)
                    elseif (aux(i+1) .le. x .and. x .le. aux(i+2)) then
                        soma = soma + VetorSol(i)*Horner(phi(i,2,:),x,0)
                    elseif (aux(i+2) .le. x .and. x .le. aux(i+3)) then
                        soma = soma + VetorSol(i)*Horner(phi(i,3,:),x,0)
                    elseif (aux(i+3) .le. x .and. x .le. aux(i+4)) then
                        soma = soma + VetorSol(i)*Horner(phi(i,4,:),x,0)
                    end if
                end if
            end do
            u = soma
        elseif (Ordem .eq. 2) then
            soma = 0.d0
            do i = 1, Ordem+PtosInt-2
                if (aux(i) .le. x .and. x .le. aux(i+Ordem)) then
                    if (aux(i) .le. x .and. x .le. aux(i+1)) then
                        soma = soma + VetorSol(i)*Horner(phi(i,1,:),x,0)
                    elseif (aux(i+1) .le. x .and. x .le. aux(i+2)) then
                        soma = soma + VetorSol(i)*Horner(phi(i,2,:),x,0)
                    end if
                end if
            end do
            u = soma
        end if
    end function u
    
    !================================================================================================
    ! AVALIACAO DA FUNCAO f
    !================================================================================================
    function f(x)
    
    real(kind=8) :: f, x
    
    if ((a .eq. 0.d0).and.(b .eq. 1.d0)) then
         f = (pi**2)*(Sin(pi*x)-9.d0*Sin(3.d0*pi*x))!Teste 1
!         f = (pi**2)*(Cos((pi/4.d0)*x))/16.d0 ! Teste 2
    else
        f = (pi**2)*(Sin(pi*x)-9.d0*Sin(3.d0*pi*x)) + (b-a)*Derk(x)-q(x)*(a+(b-a)*x)
    end if

    end function f
    
    !================================================================================================
    ! AVALIACAO DA FUNCAO k
    !================================================================================================
    function k(x)
    
    real(kind=8) :: k, x
    k = 1.d0!Teste 1
!         k = -1.d0 !Teste 2
    
    end function k
    
    !================================================================================================
    ! AVALIACAO DA DERIVADA DA FUNCAO k
    !================================================================================================
    function Derk(x)
    
    real(kind=8) :: Derk, x
    Derk = 0.d0
    
    end function Derk
    !================================================================================================
    ! AVALIACAO DA FUNCAO q
    !================================================================================================
    function q(x)
    
    real(kind=8) :: q, x
    q = 0.d0!Teste1
!     q = (pi**2)/4.d0 !Teste 2

    end function q
    
    !================================================================================================
    ! AVALIACAO DA EXATA
    !================================================================================================
    function Exata(x)
    
    real(kind=8) :: Exata, x
    
    if ((a .eq. 0.d0).and.(b .eq. 1.d0)) then
        Exata = Sin(pi*x) - Sin(3.d0*pi*x)!Teste1
!     Exata = (-1.d0/3.d0)*Cos((pi*x)/2.d0) - (sqrt(2.d0)/6.d0)*Sin((pi*x)/2.d0)+(1.d0/3.d0)*Cos((pi*x)/4.d0) !Teste 2
    else
        Exata = Sin(pi*x) - Sin(3.d0*pi*x) + a + (b-a)*x
    end if
    end function Exata
    
    !================================================================================================
    ! FUNCAO DO PRODUTO INTERNO L
    !================================================================================================
    function Funcpi(u,v,x)
    
    real(kind=8) :: Funcpi, u(:), v(:), x
    
    Funcpi = k(x)*Horner(u,x,1)*Horner(v,x,1) + q(x)*Horner(u,x,0)*Horner(v,x,0)
    
    end function Funcpi
    
    !================================================================================================
    ! FUNCAO DO PRODUTO INTERNO USUAL
    !================================================================================================
    function Funcpi_b(u,x)
    
    real(kind=8) :: Funcpi_b, u(:), x
    
    Funcpi_b = f(x)*Horner(u,x,0)
   
    end function Funcpi_b
    
    !================================================================================================
    ! PRODUTO INTERNO L
    !================================================================================================
    function ProdIntL(phii,phij,indi,indj,m,k)
        implicit none
        
        integer :: i,j,indi,indj,dif,m,k
        real(kind=8) :: phii(:,:), phij(:,:), soma, ProdIntL, inf, sup
        real(kind=8), allocatable :: aux(:) 
        
        allocate(aux(2*m+k-2))
        call Particao(m,k,Part,PartExt)
        aux = PartExt(2:2*m+k-1)
        dif = abs(indi-indj)
        soma = 0.d0
    
        do i = 1, m-dif
            inf = aux(indj+i-1)
            sup = aux(indj+i)
            soma = soma + Romberg(phii(dif+i,:),phij(i,:),inf,sup)
        end do

        ProdIntL = soma
    end function ProdIntL
    
    !================================================================================================
    ! PRODUTO INTERNO USUAL
    !================================================================================================
    function ProdIntUsual(phi,ind)
        implicit none
        
        integer :: i, j, ind
        real(kind=8) :: ProdIntUsual, soma, inf, sup
        real(kind=8) :: phi(:,:)
        real(kind=8), allocatable :: T(:), aux(:)
        
        allocate(aux(PtosInt+2))
        call Particao(Ordem,PtosInt,Part,PartExt)
        aux = PartExt(2:2*Ordem+PtosInt-1)
        soma = 0.d0
        
        do i = 1, Ordem
            inf = aux(ind+i-1)
            sup = aux(ind+i)
            soma = soma + Romberg_b(phi(i,:),inf,sup)
        end do
        
        ProdIntUsual = soma
    end function ProdIntUsual
    
    !================================================================================================
    ! INTEGRACAO DE ROMBERG PARA O PRODUTO INTERNO USUAL
    !================================================================================================
     function Romberg_b(phi,a,b)
        implicit none
        
        integer :: iter = 20, np, k, i, j
        real(kind=8) :: a, b, h, tol, m, Romberg_b
        real(kind=8) :: phi(:)
        real(kind=8), allocatable :: T(:)
        
        allocate(T(iter))
        T = 0.d0; k = 2; np = 2
        tol = 1e-5; h = b-a
        T(1) = h*((Funcpi_b(phi,a)+Funcpi_b(phi,b))/2.d0)
        T(2) = T(1)/2.d0
        h = h/2.d0
        T(2) = T(2) + h*Funcpi_b(phi,a+h)
        T(1) = (4.d0*T(2)-T(1))/3.d0
        
        do while (k .lt. iter .and. abs(T(2)-T(1)) .ge. tol)
            k = k + 1
            h = h/2.d0
            T(k) = T(k-1)/2.d0

            do i = 1, np
                T(k) = T(k) + h*Funcpi_b(phi,a+(2*i-1)*h)
            end do
            m = 1.d0
            do j = 1, k-1
                m = m*4.d0
                T(k-j) = (m*T(k-j+1)-T(k-j))/(m-1.d0)
            end do
            np = np*2

        end do
        Romberg_b = T(1)
    
     end function Romberg_b
    
    !================================================================================================
    ! INTEGRACAO DE ROMBERG PARA O PRODUTO INTERNO L
    !================================================================================================
     function Romberg(phii,phij,a,b)
        implicit none
        
        integer :: iter = 20, np, k, i, j
        real(kind=8) :: a, b, h, tol, m, Romberg
        real(kind=8) :: phii(:), phij(:)
        real(kind=8), allocatable :: T(:)
        
        allocate(T(iter))
        T = 0.d0; k = 2; np = 2
        tol = 1e-5; h = b-a
        T(1) = h*((Funcpi(phii,phij,a)+Funcpi(phii,phij,b))/2.d0)
        T(2) = T(1)/2.d0
        h = h/2.d0
        T(2) = T(2) + h*Funcpi(phii,phij,a+h)
        T(1) = (4.d0*T(2)-T(1))/3.d0
        
        do while (k .lt. iter .and. abs(T(2)-T(1)) .ge. tol)
            k = k + 1
            h = h/2.d0
            T(k) = T(k-1)/2.d0

            do i = 1, np
                T(k) = T(k) + h*Funcpi(phii,phij,a+(2*i-1)*h)
            end do
            m = 1.d0
            do j = 1, k-1
                m = m*4.d0
                T(k-j) = (m*T(k-j+1)-T(k-j))/(m-1.d0)
            end do
            np = np*2

        end do
        Romberg = T(1)
    
     end function Romberg
    
    !================================================================================================
    ! PARTICAO ESTENDIDA
    !================================================================================================
    subroutine Particao(m,k,par,pes)
        implicit none
        
        real(kind=8) :: h
        integer :: m,k,i
        real(kind=8), allocatable :: par(:), pes(:)
        
        h = 1./(k+1)
        
        do i = 1,k+2
            par(i) = (i-1)*h
            pes(i+m-1) = par(i)
        end do
        
        do i = 1,m-1
            pes(i) = (i-m)*h
            pes(i+k+m+1) = (i*h)+1
        end do
        
    end subroutine Particao
    
    !================================================================================================
    ! B-SPLINES
    !================================================================================================
    function BSplines(m,k)
        implicit none
        
        integer :: m,k,i,j,n
        real(kind=8) :: h, Q1
        real(kind=8), allocatable :: Qi(:,:),Qiaux(:,:),pol1(:),pol2(:), BSplines(:,:,:,:) 
        
        allocate(pol1(2),pol2(2),Qi(m,m),Qiaux(m,m),BSplines(2,m+k,m,m))
        call Particao(m,k,Part,PartExt)
        
        Qi = 0.d0; Qiaux = 0.d0; h = Part(2) - Part(1); Q1 = 1/h
        
        Qi(1,1:2) = (/-PartExt(1), 1.d0/)*(Q1/(Part(3)-Part(1)))
        Qi(2,1:2) = (/PartExt(3), -1.d0/)*(Q1/(Part(3)-Part(1)))
        Qiaux(1,1:2) = Translacao(Qi(1,:),h)
        Qiaux(2,1:2) = Translacao(Qi(2,:),h)
        
        BSplines(1,1,1,:) = Qi(1,:)
        BSplines(1,1,2,:) = Qi(2,:)
        
        do i = 2, m+k
            do j = 1, 2
                BSplines(1,i,j,:) = Translacao(Qi(j,:),(i-1)*h)
            end do
        end do
        
        if (m .gt. 2) then
            BSplines = 0.d0
            pol1 = (/-PartExt(1), 1.d0/)
            
            do j = 3, m
                pol2 = (/PartExt(1+j), -1.d0/)
                Qi(1,1:j) = ProdPoli(pol1,Qi(1,:))/(h*j)
                
                do i = 2, j-1
                    Qi(i,1:j) = (SomPoli(ProdPoli(Qi(i,:),pol1),ProdPoli(Qiaux(i-1,:),pol2)))/(h*j)
                end do
                Qi(j,1:j) = ProdPoli(Qiaux(j-1,:),pol2)/(h*j)
                
                do i = 1, j
                    Qiaux(i,1:j) = Translacao(Qi(i,:),h)
                end do
                
                if (j .ge. m-1) then
                    BSplines(j-m+2,1,:,:) = Qi
                    
                    do i = 2, m+k
                        do l = 1, m
                            BSplines(j-m+2,i,l,:) = Translacao(Qi(l,:),(i-1)*h)
                        end do
                    end do
                end if
            end do
        end if
        
    end function BSplines 
    
    !================================================================================================
    ! TRANSLACAO DE POLINOMIOS
    !================================================================================================
    function Translacao(p,a)
        implicit none
        
        real(kind=8) :: a
        integer :: i,n
        real(kind=8) :: p(:)
        real(kind=8), allocatable :: Translacao(:),aux(:,:),pol(:),aux1(:)
        
        n = size(p)
        allocate(aux(n,n),pol(2),Translacao(n),aux1(n))
        pol = (/-a,1.d0/); aux(1,1) = 1.d0; Translacao = 0.d0
        
        do i = 2, n
            aux(i,:) = ProdPoli(aux(i-1,:),pol)
        end do

        do i = 1, n
            aux1 = p(i)*aux(i,:)
            Translacao = SomPoli(Translacao,aux1)
        end do
        
    end function Translacao
    
    !================================================================================================
    ! SOMA DE POLINOMIOS
    !================================================================================================
    function SomPoli(p,q)
        implicit none
        
        integer :: m,n,i
        real(kind=8) :: p(:), q(:)
        real(kind=8), allocatable :: SomPoli(:)
        
        n = size(p); m = size(q)
        allocate(SomPoli(max(m,n)))
        SomPoli = 0
        
        do i = 1,min(m,n)
            SomPoli(i) = p(i)+q(i)
        end do
        
    end function SomPoli
    
    !================================================================================================
    ! PRODUTO DE POLINOMIOS
    !================================================================================================
    function ProdPoli(p,q)
        implicit none
        
        real(kind=8) :: soma, p(:),q(:)
        integer :: k, jmax, jmin, j, m, n, i
        real(kind=8), allocatable :: ProdPoli1(:), ProdPoli(:)
        
        n = size(p); m = size(q)
        allocate(ProdPoli1(m+n-1))
        
        do k = 1,m+n-1
            soma = 0
            jmax = max(1,k + 1 - m)
            jmin = min(k,n)
            
            do j = jmax,jmin
                soma = soma + p(j)*q(k - j + 1)
            end do
            
            ProdPoli1(k) = soma
        end do
        
        i = size(ProdPoli1)
        
        do while (ProdPoli1(i) .eq. 0)
            i = i - 1
        end do
        allocate(ProdPoli(i))
        ProdPoli = ProdPoli1(1:i)
        
    end function ProdPoli
    
    !================================================================================================
    ! ALGORITMO DE HORNER
    !================================================================================================
    function Horner(p,a,derivada)
        implicit none
        
        integer :: i,j,derivada,n
        real(kind=8) :: a, Horner
        real(kind=8) :: p(:)
        real(kind=8), allocatable :: aux(:)
        
        n = size(p)
        allocate(aux(n))
        aux = p

        if (derivada .le. n-1) then
            do i = 1, derivada+1
                do j = n-1, i, -1
                    aux(j) = aux(j+1)*a + aux(j)
                end do
            end do
            Horner = aux(derivada+1)*Factorial(derivada)
            else 
                Horner = 0.d0
        end if
    
    end function Horner
    
    !================================================================================================
    ! FATORIAL DE UM NUMERO
    !================================================================================================
    function Factorial(n)
        implicit none
        
        integer :: Factorial,i,n 
        Factorial = 1 
        do i= 1, n 
            Factorial = Factorial*i 
        end do
    end function Factorial
    
    !================================================================================================
    ! NORMA INFINITO
    !================================================================================================
    function Norma(u,v)
        implicit none
        integer :: n
        real(kind=8) :: u(:), v(:), Norma, aux1, aux2
        
        n = size(u); aux1 = abs(u(1)-v(1))
        do i = 2, n
            aux2 = abs(u(i)-v(i))
            Norma = max(aux1,aux2)
            aux1 = Norma
        end do
    
    end function Norma
    !================================================================================================
    ! ELIMINACAO GAUSSIANA MATRIZ DE BANDA
    !================================================================================================
    function GaussBanda(A,b,bn)
        implicit none
        integer :: n, bn, n1,n2,aux3,aux4,aux5
        integer, allocatable :: aux2(:)
        real(kind=8) :: m, A(:,:), b(:), soma
        real(kind=8), allocatable :: GaussBanda(:), Aux(:,:)
        integer :: i, j, k
        
        n = size(A(1,:)); n1 = int(bn/2)+1; n2 = int(bn/2)
        allocate(Aux(bn,n), GaussBanda(n),aux2(n2))
        Aux = 0.d0
        
        do i = 1, n1
            do j = 1, n-n1+i
                Aux(i,j) = A(j,n1-i+j)
            end do
        end do
        
        do i = n1+1, bn
            Aux(i,:) = Aux(bn-i+1,:)
        end do
    
        do i = 2, n
            do k = 1, n2
                aux2(k) = i
            end do
            aux4 = n2-1
            aux5 = 0
            
            do j = n2+2,bn
                m = Aux(j,i-1) / Aux(n2+1,i-1)
                aux3 = n2
                
                do k = j-1, j-n2,-1
                    Aux(k,aux2(aux3)) = Aux(k,aux2(aux3))-m*Aux(aux3,i-1)
                    aux3 = aux3-1
                end do
                
                do k = 1, aux4
                    aux2(k) = aux2(k)+1
                end do
                
                aux4 = aux4-1
                b(i+aux5) = b(i+aux5) - m*b(i-1)
                aux5 = aux5+1
                Aux(j,i-1) = 0.d0
            end do
        end do
        
        GaussBanda(n) = b(n) / Aux(n1,n)
        
        do i = n-1,1,-1
            soma = 0.d0
            do j = 1,n2
                soma = soma + Aux(j,i)*GaussBanda(i+n2-j+1)
            end do
            GaussBanda(i) = (b(i)-soma) / Aux(n1,i)
        end do

    end function GaussBanda

end program Projeto
