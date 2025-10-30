subroutine make_eAmatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, l, m, n
    
    !初期化
    BmatXmat = 0.0d0
    DmatBmatXmat = 0.0d0
    eAmatrix = 0.0d0
    
    !BmatXmatの作成
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dir
                do j = 1, dir
                    do m = 1, scount*dim
                        BmatXmat(k, l, i, j) = BmatXmat(k, l, i, j) + Bmatrix(k, l, i, m) * eXmatrix(k, m, j)
                    end do
                end do
            end do
        end do
    end do
    
    !DmatBmatXmatの作成
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dir
                do j = 1, dir
                    do m = 1, dir
                        DmatBmatXmat(k, l, i, j) = DmatBmatXmat(k, l, i, j) + Dmatrix(k, i, m) * BmatXmat(k, l, m, j)
                    end do
                end do
            end do
        end do
    end do
    
    !eAmatrixの作成(先生博論式(2.15))
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dir
                do j = 1, dir
                    eAmatrix(k, l, i, j) = eAmatrix(k, l, i, j) + DmatBmatXmat(k, l, i, j) + Dmatrix(k, i, j)
                end do
            end do
        end do
    end do

    !eAmatrixの出力
    open(1, file = 'eAmatrix.dat', status = 'replace')
        do k = 1, ecount
            do l = 1, ipn
                write(1, '(A, I2, A, I2)') "element_", k, "-", l
                do i = 1, 6
                    write(1, '(6F15.7)') (eAmatrix(k, l, i, j), j = 1, 6)
                end do
            end do
        end do
    close(1)
    
    print *, 'made eAmatrix(output file : "eAmatrix.dat")'


end subroutine