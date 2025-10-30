subroutine make_Fmatrix

    use set_parameter
    
    implicit none
    
    integer i, i1, j, j1, k
    
    call make_eFmatrix
    
    !初期化
    Fmatrix =0.0d0
    
    !作成
    do k = 1, ecount
        do j = 1, 6
            do i1 = 1, ipn          !積分点番号の行
                do j1 = 1, dim      !積分点番号の次元
                    Fmatrix((element(k, i1) - 1) * dim + j1, j) = Fmatrix((element(k, i1) - 1) * dim + j1, j) + eFmatrix(k, (i1 - 1) * dim + j1, j)
                end do
            end do
        end do
    end do
    
    Fmatrix = -Fmatrix
    
    !Fマトリックスの出力
    open(1, file = 'Fmatrix.dat', status = 'replace')
        do i=1, ncount*dim
            write(1, '(6F15.7)') (Fmatrix(i, j), j = 1, 6)
        end do
    close(1)

    print *, 'made Fmatrix (outputfile : "Fmatrix.dat")'

end subroutine