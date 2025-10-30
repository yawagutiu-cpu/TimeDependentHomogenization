subroutine solve_equation(A, B, solution, m, n, nrhs)    !LAPACKで連立方程式[A][X]=[B]の解を計算する関数
                                                            !入力値([A], [B])を変化させることなく解[X]を返す 
                                                            !入力([A], [B], [X], ubound([A], 1), ubound([A], 2), ubound([B], 2))
                                                            !出力([X])
                                       
    use set_parameter
    
    implicit none
    
    integer i, j, lda, ldb
    integer, intent(in) :: m, n, nrhs
    double precision, intent(in) :: A(m, n), B(m, nrhs)   !引数の指定(メモリを指定する必要あり)
    double precision, intent(inout) :: solution(m, nrhs)
    double precision A_temp(m, n)
    
    !LAPACKでの計算準備
    lda = m
    ldb = ubound(B, 1)
    allocate (ipiv(m))
    ipiv = 0
    
    !入力値([A], [B])を保存
    A_temp = 0.0d0
    
    do i = 1, m
        do j = 1, n
            A_temp(i, j) = A(i, j)
        end do
    end do
    
    do i = 1, ldb
        do j = 1, nrhs
            solution(i, j) = B(i, j)
        end do
    end do
    
    !求解
    call dgetrf(m, n, A_temp, lda, ipiv, info)
    call dgetrs(trans, n, nrhs, A_temp, lda, ipiv, solution, ldb, info)
    
    !配列サイズクリア
    deallocate (ipiv)
    
end subroutine