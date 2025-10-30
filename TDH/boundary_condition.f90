subroutine boundary_condition

    use set_parameter
    
    implicit none
    
    integer i, j, k, m, n, lda
    
    !節点1(原点)のx,y,z方向を固定
    do i = 1, dim
        Kmatrix(i, i) = Kmatrix(i, i) + lambda
    end do
    
    do i = 1, pair
        do j = 1, dim
        Kmatrix(dim*penalty_number(i, 1)+j-dim, dim*penalty_number(i, 1)+j-dim) = Kmatrix(dim*penalty_number(i, 1)+j-dim, dim*penalty_number(i, 1)+j-dim) + lambda
        Kmatrix(dim*penalty_number(i, 1)+j-dim, dim*penalty_number(i, 2)+j-dim) = Kmatrix(dim*penalty_number(i, 1)+j-dim, dim*penalty_number(i, 2)+j-dim) - lambda
        Kmatrix(dim*penalty_number(i, 2)+j-dim, dim*penalty_number(i, 1)+j-dim) = Kmatrix(dim*penalty_number(i, 2)+j-dim, dim*penalty_number(i, 1)+j-dim) - lambda
        Kmatrix(dim*penalty_number(i, 2)+j-dim, dim*penalty_number(i, 2)+j-dim) = Kmatrix(dim*penalty_number(i, 2)+j-dim, dim*penalty_number(i, 2)+j-dim) + lambda
        end do
    end do
    
    !ペナルティ法を適用した後のKmatrixを出力
    open(1, file = 'Kmatrix_penalty.dat', status = 'replace')
        do i = 1, ncount*dim
            write(1, '(192F25.7)') (Kmatrix(i, j), j = 1, ncount*dim)
        end do
    close(1)
    
    !LU分解しておく（時間依存の関数で使う）
    m = ubound(Kmatrix, 1)
    n = ubound(Kmatrix, 2)
    lda = m
    allocate (ipiv_k(m))
    ipiv_k = 0
     
    Kmat_LU = Kmatrix
     
    call dgetrf(m, n, Kmat_LU, lda, ipiv_k, info)
    
end subroutine