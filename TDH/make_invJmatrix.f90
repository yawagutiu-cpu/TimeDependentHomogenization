subroutine make_invJmatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, l, m, n, nrhs
    
    !初期化
    Jmatrix = 0.0d0
    invJmatrix = 0.0d0
    
    !Jmatrixの作成
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dim
                do j = 1, dim
                    do m = 1, ipn
                        Jmatrix(k, l, i, j) = Jmatrix(k, l, i, j) + Hmatrix(k, l, i, m) * node(element(k, m), j)
                    end do
                end do
            end do
        end do
    end do
    
    !det_Jの計算
    det_J = 0.0d0
    
     do k = 1, ecount
        do l = 1, ipn
            det_J(k, l) = Jmatrix(k, l, 1, 1) * Jmatrix(k, l, 2, 2) * Jmatrix(k, l, 3, 3) + Jmatrix(k, l, 1, 2) * Jmatrix(k, l, 2, 3) * Jmatrix(k, l, 3, 1) + Jmatrix(k, l, 1, 3) * Jmatrix(k, l, 2, 1) * Jmatrix(k, l, 3, 2) - Jmatrix(k, l, 1, 3) * Jmatrix(k, l, 2, 2) * Jmatrix(k, l, 3, 1) - Jmatrix(k, l, 1, 2) * Jmatrix(k, l, 2, 1) * Jmatrix(k, l, 3, 3) - Jmatrix(k, l, 1, 1) * Jmatrix(k, l, 3, 2) * Jmatrix(k, l, 2, 3)
        end do
     end do
   
    !invJmatrixの作成
    !配列サイズの指定
    allocate (Imatrix(dim, dim))
    allocate (temp(dim, dim))
    allocate (invtemp(dim, dim))
    
    !初期化
    Imatrix = 0.0d0
    temp = 0.0d0
    invtemp = 0.0d0
    
    !単位行列の作成    
    do i = 1, dim
        Imatrix(i, i) = 1.0d0
    end do
    
    !計算(各要素，各積分点でinvJmatrixを求める)
    do k = 1, ecount
        do l = 1, ipn
        
            !k要素l積分点のJmatrixをtempに一時的に保存
            do i = 1, dim
                do j = 1, dim
                    temp(i, j) = Jmatrix(k, l, i, j)
                end do
            end do
            
            ![J][X]=[E]を解く
            m = ubound(temp, 1)
            n = ubound(temp, 2)
            nrhs = ubound(Imatrix, 2)
            call solve_equation(temp, Imatrix, invtemp, m, n, nrhs)    !入力([A], [B], [X], ubound([A], 1), ubound([A], 2), ubound([B], 2))   
                                                                                        !tempとImatrixは変更されない
            
            !求めたinvtempをk要素l積分点のinvJmatrixとして保存
            do i = 1, dim
                do j = 1, dim
                    invJmatrix(k, l, i, j) = invtemp(i, j)
                end do
            end do
            
        end do
    end do
    
    !配列サイズクリア
    deallocate (Imatrix)
    deallocate (temp)
    deallocate (invtemp)

   
   !出力
   open(1, file = 'Jmatrix.dat', status = 'replace')
        do k = 1, ecount
            do l = 1, ipn
                write(1, '(A, I2, A, I2)') "element_", k, "-", l
                do i = 1, dim
                    write(1, '(3F15.7)') (Jmatrix(k, l, i, j), j = 1, dim)
                end do
            end do
        end do
    close(1)
    
    print *, 'made Jmatrix(output file : "Jmatrix.dat")'
    
    open(2, file = 'invJmatrix.dat', status = 'replace')
        do k = 1, ecount
            do l = 1, ipn
                write(2, '(A, I2, A, I2)') "element_", k, "-", l
                do i = 1, dim
                    write(2, '(3F15.7)') (invJmatrix(k, l, i, j), j = 1, dim)
                end do
            end do
        end do
    close(2)
    
    print *, 'made invJmatrix(output file : "invJmatrix.dat")'

end subroutine