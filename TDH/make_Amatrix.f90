subroutine make_Amatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, l, m, n, nrhs
    
    !初期化
    Amatrix = 0.0d0
    invAmatrix = 0.0d0

    !Amatrixの作成(正確にはAマトリックスの体積平均(巨視的弾性剛性テンソル))
    call make_eAmatrix
    
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dir
                do j = 1, dir
                    Amatrix(i, j) = Amatrix(i, j) + eAmatrix(k, l, i, j) * det_J(k, l)
                end do
            end do
        end do
    end do
      
    !均質化
    Amatrix = Amatrix / Vol
    
    !invAmatrixを求める
    !単位行列の用意
    allocate (Imatrix(dir, dir))
    
    Imatrix = 0.0d0
    
    do i = 1, dir
        Imatrix(i, i) = 1.0d0
    end do
    
    !計算
    m = ubound(Amatrix, 1)
    n = ubound(Amatrix, 2)
    nrhs = ubound(Imatrix, 2)
    call solve_equation(Amatrix, Imatrix, invAmatrix, m, n, nrhs)
    
    !配列サイズクリア
    deallocate (Imatrix)
    
    !出力
    open(1, file = 'Amatrix.dat', status = 'replace')
        do i = 1, 6
            write(1, '(6F15.7)') (Amatrix(i, j), j = 1, 6)
        end do
    close(1)
    
    print *, 'made Amatrix(output file : "Amatrix.dat")'
    
    open(2, file = 'invAmatrix.dat', status = 'replace')
        do i = 1, 6
            write(2, '(6F15.7)') (invAmatrix(i, j), j = 1, 6)
        end do
    close(2)
    
    print *, 'made invAmatrix(output file : "invAmatrix.dat")'
    
    !均質化された弾性定数の出力
    open(3, file = 'Effective_Modules.dat', status = 'replace')
        write(3, '(A)') "symmetry : 直交異方性"
        write(3, '(A, F15.7)'), "E11 :", 1.0d0 / invAmatrix(1, 1)
        write(3, '(A, F15.7)'), "E22 :", 1.0d0 / invAmatrix(2, 2)
        write(3, '(A, F15.7)'), "E33 :", 1.0d0 / invAmatrix(3, 3)
        write(3, '(A, F15.7)'), "ν12 :", -invAmatrix(1, 2) / invAmatrix(1, 1)
        write(3, '(A, F15.7)'), "ν13 :", -invAmatrix(1, 3) / invAmatrix(1, 1)
        write(3, '(A, F15.7)'), "ν21 :", -invAmatrix(2, 1) / invAmatrix(2, 2)
        write(3, '(A, F15.7)'), "ν23 :", -invAmatrix(2, 3) / invAmatrix(2, 2)
        write(3, '(A, F15.7)'), "ν31 :", -invAmatrix(3, 1) / invAmatrix(3, 3)  
        write(3, '(A, F15.7)'), "ν32 :", -invAmatrix(3, 2) / invAmatrix(3, 3)  
        write(3, '(A, F15.7)'), "G12 :", 1.0d0 / invAmatrix(4, 4)  
        write(3, '(A, F15.7)'), "G23 :", 1.0d0 / invAmatrix(5, 5)  
        write(3, '(A, F15.7)'), "G31:", 1.0d0 / invAmatrix(6, 6)
    close(3)
    
end subroutine