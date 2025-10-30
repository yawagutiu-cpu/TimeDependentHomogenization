subroutine make_eFmatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, l, m
    
    !初期化
    eFmatrix = 0.0d0
    eFmat_c = 0.0d0
    eFmat_n = 0.0d0
    eKcn_inveKnn_eFn = 0.0d0
    
    ![eF] = [B]t * [D] * |J| * w (w=1)
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, scount*dim
                do j = 1, 6
                    eFmatrix(k, i, j) = eFmatrix(k, i, j) + tBmatDmat(k, l, i, j) * det_J(k, l)
                end do
            end do 
        end do
    end do
    
    !非適合要素を使用する場合の処理*****************************************************************************************
    if(mode == 2) then
        !関連しているモード（適合か非適合）ごとに行列を下記のように分ける
        !------------------------------
        !              |              |
        !              |     eFc      |
        !              |    (24*6)    |
        ! eFmatrix =   |              |
        !              |--------------|
        !              |     eFn      |
        !              |    (9*6)     |
        !------------------------------    
        do k = 1, ecount
            !eFc
            do i = 1, ipn*dim
                do j = 1, dir
                    eFmat_c(k, i, j) = eFmatrix(k, i, j)
                end do
            end do
            !eFn
            do i = 1, dim*dim
                do j = 1, dir
                    eFmat_n(k, i, j) = eFmatrix(k, ipn*dim + i, j)
                end do
            end do
        end do
        
        !eFc, eFnの出力
        open(1, file = 'eFmat_c.dat', status = 'replace')
            do k = 1, ecount
                write(1, '("element", I2)') k
                do i = 1, ipn*dim
                    write(1, '(6F15.7)') (eFmat_c(k, i, j), j = 1, dir)
                end do
            end do
        close(1)
        
        open(2, file = 'eFmat_n.dat', status = 'replace')
            do k = 1, ecount
                write(2, '("element", I2)') k
                do i = 1, dim*dim
                    write(2, '(6F15.7)') (eFmat_n(k, i, j), j = 1, dir)
                end do
            end do
        close(2)
        
        print *, 'distributed eFmatrix'
        
        !eKmatrixの時と同様に非適合モードに対応する変位場（u9,u10,u11,v9,v10,v11,w9,w10,w11）に関する式を消去する
        ![eKcn]*[eKnn]^-1*[eFn]の計算
        do k = 1, ecount
            do i = 1, ipn*dim
                do j = 1, dir
                    do m = 1, dim*dim
                        eKcn_inveKnn_eFn(k, i, j) = eKcn_inveKnn_eFn(k, i, j) + eKcn_inveKnn(k, i, m) * eFmat_n(k, m, j)
                    end do
                end do
            end do
        end do
        
        ![eF]の作成
        !便宜上eFmatrixに上書きする
        eFmatrix = 0.0d0
        
        do k = 1, ecount
            do i = 1, ipn*dim
                do j = 1, dir
                    eFmatrix(k, i, j) = eFmat_c(k, i, j) - eKcn_inveKnn_eFn(k, i, j)
                end do
            end do
        end do
    end if
    !********************************************************************************************************************
    
    !eKmatrixの出力
    open(1, file = 'eFmatrix.dat', status = 'replace')
        do k = 1, ecount
            write(1, '(A, I2)') "element_", k
            do i = 1, scount*dim
                write(1, '(6F15.7)') (eFmatrix(k, i, j), j = 1, 6)
            end do
        end do
    close(1)
        
    print *, 'made eFmatrix (output file : "eFmatrix.dat")'

end subroutine