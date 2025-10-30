subroutine make_Hmatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, l
    
    !初期化
    Hmatrix = 0.0d0
    
    !Hmatrixの作成（形状関数：Ni = (1+ξ*ξi)*(1+η*ηi)*(1+ζ*ζi)/8の偏微分）
    do k = 1, ecount
        do i = 1, ipn
            do j = 1, ipn
            !適合要素形状関数
            !polarで形状関数（i=1～8）を決定し，gaussを代入することでガウス・ルシャンドル積分に使えるようにする
                !ξで偏微分
                Hmatrix(k, j, 1, i) = polar(i, 1) * (1.0d0 + gauss(j, 2) * polar(i, 2)) * (1.0d0 + gauss(j, 3) * polar(i, 3)) / 8.0d0
                !ηで偏微分
                Hmatrix(k, j, 2, i) = (1.0d0 + gauss(j, 1) * polar(i, 1)) * polar(i, 2) * (1.0d0 + gauss(j, 3) * polar(i, 3)) / 8.0d0
                !ζで偏微分
                Hmatrix(k, j, 3, i) = (1.0d0 + gauss(j, 1) * polar(i, 1)) * (1.0d0 + gauss(j, 2) * polar(i, 2)) * polar(i, 3) / 8.0d0
            end do
            if(mode == 2) then
            !非適合要素形状関数
                !ξで微分
                Hmatrix(k, i, 1, 9) = -2.0d0 * gauss(i, 1)
                !ηで微分
                Hmatrix(k, i, 2, 10) = -2.0d0 * gauss(i, 2)
                !ζで微分
                Hmatrix(k, i, 3, 11) = -2.0d0 * gauss(i, 3)
            end if
        end do
    end do
    
    !出力
    open(1, file = 'Hmatrix.dat', status = 'replace')
        do k = 1, ecount
            do l = 1, ipn
                write(1, '(A, I2, A, I2)') "element_", k, "-", l
                do i = 1, dim
                    if(mode == 1) then
                        write(1, '(8F15.7)') (Hmatrix(k, l, i, j), j = 1, scount)
                    else if(mode == 2) then
                        write(1, '(11F15.7)') (Hmatrix(k, l, i, j), j = 1, scount)
                    end if
                end do
            end do
        end do
    close(1)
    
    print *, 'made Hmatrix(output file : "Hmatrix.dat")'

end subroutine