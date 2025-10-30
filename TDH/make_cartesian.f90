subroutine make_cartesian

    use set_parameter 
    
    implicit none
    
    integer i, j, k, l, m
    
    !初期化
    cartesian = 0.0d0
    
    !作成
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dim
                do j = 1, scount    
                    do m = 1, dim
                        cartesian(k, l, i, j) = cartesian(k, l, i, j) + invJmatrix(k, l, i, m) * Hmatrix(k, l, m, j)
                    end do
                end do
            end do
        end do
    end do
    
    !出力
     open(1, file = 'cartesian.dat', status = 'replace')
        do k = 1, ecount
            do l = 1, ipn
                write(1, '(A, I2, A, I2)') "element_", k, "-", l
                do i = 1, dim
                    if(mode == 1) then
                        write(1, '(8F15.7)') (cartesian(k, l, i, j), j = 1, scount)
                    else if(mode == 2) then
                        write(1, '(11F15.7)') (cartesian(k, l, i, j), j = 1, scount)
                    end if
                end do
            end do
        end do
    close(1)
    
    print *, 'made cartesian(output file : "cartesian.dat")'

end subroutine