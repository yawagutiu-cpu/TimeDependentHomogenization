subroutine make_Xmatrix

    use set_parameter

    implicit none

    ! 変数宣言
    integer :: i, j, k, l, nrhs, lda, ldb
    integer :: m, n

    double precision, allocatable :: A_temp(:,:)


    ! 初期化
    Xmatrix = 0.0d0
    eXmat_c = 0.0d0
    eKnc_eXc = 0.0d0
    eXmat_n = 0.0d0
    eXmatrix = 0.0d0

    ! 必要な行列の作成
    call make_Fmatrix

    ! 行列サイズの取得
    m = ubound(Kmatrix, 1)
    n = ubound(Kmatrix, 2)
    nrhs = ubound(Fmatrix, 2)
    lda = m
    ldb = m

    ! 配列の動的割り当て
    allocate(ipiv(m))
    allocate(A_temp(m, n))

    ! デバッグ用出力
    print *, 'OK'
    print *, 'Dimensions of Kmatrix:', m, n
    print *, 'Dimensions of Fmatrix:', nrhs

    ! Kmatrix をコピーして A_temp に保持
    A_temp = Kmatrix
    Xmatrix = Fmatrix

    ! 求解準備（LU分解 & 前進後退代入）
    call dgetrf(m, n, A_temp, lda, ipiv, info)
    if (info /= 0) then
        print *, 'Error in dgetrf: info =', info
        stop
    end if

    call dgetrs('N', n, nrhs, A_temp, lda, ipiv, Xmatrix, ldb, info)
    if (info /= 0) then
        print *, 'Error in dgetrs: info =', info
        stop
    end if

    print *, 'Solution computed successfully'

    ! メモリ解放
    deallocate(ipiv, A_temp)

    ! Xmatrix の出力
    open(1, file = 'Xmatrix.dat', status = 'replace')
    do i = 1, ncount * dim
        write(1, '(6F15.7)') (Xmatrix(i, j), j = 1, dir)
    end do
    close(1)
    print *, 'made Xmatrix (output file："Xmatrix.dat")'

    ! 適合要素を使用する場合
    if (mode == 1) then
        ! χmatrix を eχmatrix に再配置
        do k = 1, ecount
            do l = 1, ipn
                do i = 1, dim
                    do j = 1, dir
                        eXmatrix(k, dim * (l - 1) + i, j) = Xmatrix(dim * (element(k, l) - 1) + i, j)
                    end do
                end do
            end do
        end do

    ! 非適合要素を使用する場合
    else if (mode == 2) then
        ! Xmatrix を要素ごとに分配して eXmat_c を生成
        do k = 1, ecount
            do l = 1, ipn
                do i = 1, dim
                    do j = 1, dir
                        eXmat_c(k, (l - 1) * dim + i, j) = Xmatrix((element(k, l) - 1) * dim + i, j)
                    end do
                end do
            end do
        end do

        ! eXmat_c から eXmat_n を計算
        do k = 1, ecount
            do i = 1, dim * dim
                do j = 1, dir
                    do m = 1, ipn * dim
                        eKnc_eXc(k, i, j) = eKnc_eXc(k, i, j) + eKmat_nc(k, i, m) * eXmat_c(k, m, j)
                    end do
                end do
            end do

            do i = 1, dim * dim
                do j = 1, dir
                    do m = 1, dim * dim
                        eXmat_n(k, i, j) = eXmat_n(k, i, j) + inveKmat_nn(k, i, m) * (eFmat_n(k, m, j) - eKnc_eXc(k, m, j))
                    end do
                end do
            end do
        end do

        ! eXmat_c, eXmat_n を基に eXmatrix を作成
        do k = 1, ecount
            do i = 1, ipn * dim
                do j = 1, dir
                    eXmatrix(k, i, j) = eXmat_c(k, i, j)
                end do
            end do

            do i = 1, dim * dim
                do j = 1, dir
                    eXmatrix(k, ipn * dim + i, j) = eXmat_n(k, i, j)
                end do
            end do
        end do
    end if

    ! eXmatrix の出力
    open(2, file = 'eXmatrix.dat', status = 'replace')
    do k = 1, ecount
        write(2, '(A, I2)') "element_", k
        do i = 1, scount * dim
            write(2, '(6F15.7)') (eXmatrix(k, i, j), j = 1, dir)
        end do
    end do
    close(2)
    print *, 'made eXmatrix (output file："eXmatrix.dat")'

end subroutine