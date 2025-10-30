subroutine make_Smatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, m, n, nrhs
    
    !初期化
    Smatrix = 0.0d0
    invSmatrix = 0.0d0
    
    !Smatrixの作成
    !1軸方向炭素繊維
    Smatrix(1, 1, 1) = 1.0d0 / Ef_3
    Smatrix(1, 1, 2) = -poissonf_31 / Ef_3
    Smatrix(1, 1, 3) = Smatrix(1, 1, 2)
    Smatrix(1, 2, 1) = Smatrix(1, 1, 2)
    Smatrix(1, 2, 2) = 1.0d0 / Ef_1
    Smatrix(1, 2, 3) = -poissonf_12 / Ef_1
    Smatrix(1, 3, 1) = Smatrix(1, 1, 2)
    Smatrix(1, 3, 2) = Smatrix(1, 2, 3)
    Smatrix(1, 3, 3) = Smatrix(1, 2, 2)
    
    Smatrix(1, 4, 4) = 1.0d0 / Gf_31
    Smatrix(1, 5, 5) = 2.0d0 * (Smatrix(1, 2, 2) - Smatrix(1, 2, 3))
    Smatrix(1, 6, 6) = Smatrix(1, 4, 4)
    
    !2軸方向炭素繊維
    Smatrix(2, 1, 1) = 1.0d0 / Ef_1
    Smatrix(2, 1, 2) = -poissonf_31 / Ef_3
    Smatrix(2, 1, 3) = -poissonf_12 / Ef_1
    Smatrix(2, 2, 1) = Smatrix(2, 1, 2)
    Smatrix(2, 2, 2) = 1.0d0 / Ef_3
    Smatrix(2, 2, 3) = Smatrix(2, 1, 2)
    Smatrix(2, 3, 1) = Smatrix(2, 1, 3)
    Smatrix(2, 3, 2) = Smatrix(2, 1, 2)
    Smatrix(2, 3, 3) = Smatrix(2, 1, 1)
    
    Smatrix(2, 4, 4) = 1.0d0 / Gf_31
    Smatrix(2, 5, 5) = Smatrix(2, 4, 4)
    Smatrix(2, 6, 6) = 2.0d0 * (Smatrix(2, 3, 3) - Smatrix(2, 3, 1))
    
    !3軸方向炭素繊維
    Smatrix(3, 1, 1) = 1.0d0 / Ef_1
    Smatrix(3, 1, 2) = -poissonf_12 / Ef_1
    Smatrix(3, 1, 3) = -poissonf_31 / Ef_3
    Smatrix(3, 2, 1) = Smatrix(3, 1, 2)
    Smatrix(3, 2, 2) = Smatrix(3, 1, 1)
    Smatrix(3, 2, 3) = Smatrix(3, 1, 3)
    Smatrix(3, 3, 1) = Smatrix(3, 1, 3)
    Smatrix(3, 3, 2) = Smatrix(3, 1, 3)
    Smatrix(3, 3, 3) = 1.0d0 / Ef_3
    
    Smatrix(3, 4, 4) = 2.0d0 * (Smatrix(3, 1, 1) - Smatrix(3, 1, 2))
    Smatrix(3, 5, 5) = 1.0d0 / Gf_31
    Smatrix(3, 6, 6) = Smatrix(3, 5, 5)
    
    !LAPACKを用いてinvSmatrixを求める
    !配列サイズの指定
    allocate (Imatrix(dir, dir))
    allocate (temp(dir, dir))
    allocate (invtemp(dir, dir))
    
    !初期化　
    Imatrix = 0.0d0
    temp = 0.0d0
    invtemp = 0.0d0
    
    !単位行列の用意
    do i = 1, dir
        Imatrix(i, i) = 1.0d0
    end do
    
    !各軸方向炭素繊維でそれぞれ計算する
    do k = 1, dim
        
        !k軸方向炭素繊維のSmatrixをtempに一時的に保存
        do i = 1, dir
            do j = 1, dir
                temp(i,j) = Smatrix(k, i, j)
            end do
        end do 
        
        ![S][X]=[E]を解く
        m = ubound(temp, 1)
        n = ubound(temp, 2)
        nrhs = ubound(Imatrix, 2)
        call solve_equation(temp, Imatrix, invtemp, m, n, nrhs)    !入力([A], [B], [X], ubound([A], 1), ubound([A], 2), ubound([B], 2))
                                                                    !tempとImatrixは変更されない
        
        !計算結果をk軸方向炭素繊維のinvSmatrixに保存
        do i = 1, dir
            do j = 1, dir
                invSmatrix(k, i, j) = invtemp(i, j)
            end do
        end do
  
    end do
    
    !配列サイズクリア
    deallocate (Imatrix)
    deallocate (temp)
    deallocate (invtemp)
    
end subroutine