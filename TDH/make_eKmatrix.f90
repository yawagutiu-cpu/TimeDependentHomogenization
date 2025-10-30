subroutine make_eKmatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, l, m, n, nrhs
    
    !初期化
    tBmatDmat = 0.0d0
    tBmatDmatBmat = 0.0d0
    eKmatrix = 0.0d0
    eKmat_cc = 0.0d0
    eKmat_cn = 0.0d0
    eKmat_nc = 0.0d0
    eKmat_cn = 0.0d0
    eKmat_nn = 0.0d0
    inveKmat_nn = 0.0d0
    eKcn_inveKnn = 0.0d0
    eKcn_inveKnn_eKnc = 0.0d0
    
    ![B]t[D]の作成
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, scount*dim
                do j = 1, 6
                    do m = 1, 6
                        tBmatDmat(k, l, i, j) = tBmatDmat(k, l, i, j) + Bmatrix(k, l, m, i) * Dmatrix(k, m, j)
                    end do
                end do
            end do
        end do 
    end do
        
    ![B]t[D][B]の作成
     do k = 1, ecount
        do l = 1, ipn
            do i = 1, scount*dim
                do j = 1, scount*dim
                    do m = 1, 6
                        tBmatDmatBmat(k, l, i, j) = tBmatDmatBmat(k, l, i, j) + tBmatDmat(k, l, i, m) * Bmatrix(k, l, m, j)
                    end do
                end do
            end do
        end do 
    end do
        
    !要素剛性行列の作成
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, scount*dim
                do j = 1, scount*dim
                    eKmatrix(k, i, j) = eKmatrix(k, i, j) + tBmatDmatBmat(k, l, i, j) * det_J(k, l)
                end do
            end do
        end do
    end do
    
    !非適合要素を使用する場合の処理**************************************************************************************************
    if(mode == 2) then
        !関連しているモード（適合か非適合）ごとに行列を下記のように分ける
        !--------------------------------------------
        !                |                |          |
        !                |                |   eKcn   |
        !                |  eKcc(24*24)   |          |
        ! eKmatrix  =    |                |  (24*9)  |          
        !                |                |          |
        !                |----------------|----------|
        !                |  eKnc(9*24)    |eKnn(9*9) |
        !--------------------------------------------   
        do k = 1, ecount
            do i = 1, ipn*dim
                !eKcc
                do j = 1, ipn*dim
                    eKmat_cc(k, i, j) = eKmatrix(k, i, j)
                end do
                !eKcn
                do j = 1, dim*dim
                    eKmat_cn(k, i, j) = eKmatrix(k, i, ipn*dim + j)
                end do
            end do
            
            do i = 1, dim*dim
                !eKnc
                do j = 1, ipn*dim
                    eKmat_nc(k, i, j) = eKmatrix(k, ipn*dim + i, j)
                end do
                !eKnn
                do j = 1, dim*dim
                    eKmat_nn(k, i, j) = eKmatrix(k, ipn*dim + i, ipn*dim + j)
                end do
            end do
        end do
        
        !eKcc, eKcn, eKnc, eKnnの出力
        open(1, file = 'eKmat_cc.dat', status = 'replace')
            do k = 1, ecount
                write(1, '("element_", I2)') k
                do i = 1, ipn*dim
                    write(1, '(24F15.7)') (eKmat_cc(k, i, j), j = 1, ipn*dim)
                end do
            end do
        close(1)
        
        open(2, file = 'eKmat_cn.dat', status = 'replace')
            do k = 1, ecount
                write(2, '("element_", I2)') k
                do i = 1, ipn*dim
                    write(2, '(9F15.7)') (eKmat_cn(k, i, j), j = 1, dim*dim)
                end do
            end do
        close(2)
        
        open(3, file = 'eKmat_nc.dat', status = 'replace')
            do k = 1, ecount
                write(3, '("element_", I2)') k
                do i = 1, dim*dim
                    write(3, '(24F15.7)') (eKmat_nc(k, i, j), j = 1, ipn*dim)
                end do
            end do
        close(3)
        
        open(4, file = 'eKmat_nn.dat', status = 'replace')
            do k = 1, ecount
                write(4, '("element_", I2)') k
                do i = 1, dim*dim
                    write(4, '(9F15.7)') (eKmat_nn(k, i, j), j = 1, dim*dim)
                end do
            end do
        close(4)
        
        print *, 'distributed eKmatrix'
        
        !非適合モードに対応する変位場（u9,u10,u11,v9,v10,v11,w9,w10,w11）に関する式を消去する（非適合要素テキストp252）
        !eKnnの逆行列を計算
        !配列サイズの指定
        allocate (Imatrix(dim*dim, dim*dim))
        allocate (temp(dim*dim, dim*dim))
        allocate (invtemp(dim*dim, dim*dim))
        
        !初期化
        Imatrix = 0.0d0
        temp = 0.0d0
        invtemp = 0.0d0
        
        !単位行列の作成    
        do i = 1, dim*dim
            Imatrix(i, i) = 1.0d0
        end do
        
        !計算
        do k = 1, ecount
            
            !k要素のeKmat_nnをtempに一時的に保存
            do i = 1, dim*dim
                do j = 1, dim*dim
                    temp(i, j) = eKmat_nn(k, i, j)
                end do
            end do
                
            ![eKnn][X]=[E]を解く
            m = ubound(temp, 1)
            n = ubound(temp, 2)
            nrhs = ubound(Imatrix, 2)
            call solve_equation(temp, Imatrix, invtemp, m, n, nrhs)    !入力([A], [B], [X], ubound([A], 1), ubound([A], 2), ubound([B], 2))   
                                                                        !tempとImatrixは変更されない
                
            !求めたinvtempをk要素のinveKmat_nnとして保存
            do i = 1, dim*dim
                do j = 1, dim*dim
                    inveKmat_nn(k, i, j) = invtemp(i, j)
                end do
            end do
                
        end do
        
        !配列サイズクリア
        deallocate (Imatrix)
        deallocate (temp)
        deallocate (invtemp)
                
        ![eKcn][eKnn]^-1を計算
        do k = 1, ecount
            do i = 1, ipn*dim
                do j = 1, dim*dim
                    do m = 1, dim*dim
                        eKcn_inveKnn(k, i, j) = eKcn_inveKnn(k, i, j) + eKmat_cn(k, i, m) * inveKmat_nn(k, m, j)
                    end do
                end do
            end do
        end do
        
        ![eKcn][eKnn]^-1[eKnc]を計算
        do k = 1, ecount
            do i = 1, ipn*dim
                do j = 1, ipn*dim
                    do m = 1, dim*dim
                        eKcn_inveKnn_eKnc(k, i, j) = eKcn_inveKnn_eKnc(k, i, j) + eKcn_inveKnn(k, i, m) * eKmat_nc(k, m, j)
                    end do
                end do
            end do
        end do
        
        !非適合モードに対応する変位場を消去した剛性方程式における要素剛性行列の作成
        eKmatrix = 0.0d0    !便宜上eKmatrixに上書きする
        
        do k = 1, ecount
            do i = 1, ipn*dim
                do j = 1, ipn*dim
                    eKmatrix(k, i, j) = eKmat_cc(k, i, j) - eKcn_inveKnn_eKnc(k, i, j)  !非適合要素テキスト式(104)
                end do
            end do
        end do
    end if
    !************************************************************************************************************************************
    
    !eKmatrixの出力
    open(5, file = 'eKmatrix.dat', status = 'replace')
        do k = 1, ecount
            write(5, '(A, I2)') "element_", k
            do i = 1, scount*dim
                    if(mode == 1) then
                        write(5, '(24F15.7)') (eKmatrix(k, i, j), j = 1, scount*dim)
                    else if(mode == 2) then
                        write(5, '(33F15.7)') (eKmatrix(k, i, j), j = 1, scount*dim)
                    end if
            end do
        end do
    close(5)
    
    print *, 'made eKmatrix (output file : "eKmatrix.dat")'    

end subroutine