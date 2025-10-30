subroutine make_Bmatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, l
    
    !初期化
    Bmatrix = 0.0d0
    
    !作成(cartesian行列の再配置)
    do k = 1, ecount
        do l = 1, ipn
            do j = 1, scount
                Bmatrix(k, l, 1, 3*j-2) = cartesian(k, l, 1, j)
                Bmatrix(k, l, 4, 3*j-1) = cartesian(k, l, 1, j)
                Bmatrix(k, l, 6, 3*j)   = cartesian(k, l, 1, j)
                
                Bmatrix(k, l, 4, 3*j-2) = cartesian(k, l, 2, j)
                Bmatrix(k, l, 2, 3*j-1) = cartesian(k, l, 2, j)
                Bmatrix(k, l, 5, 3*j)   = cartesian(k, l, 2, j)
                
                Bmatrix(k, l, 6, 3*j-2) = cartesian(k, l, 3, j)
                Bmatrix(k, l, 5, 3*j-1) = cartesian(k, l, 3, j)
                Bmatrix(k, l, 3, 3*j)   = cartesian(k, l, 3, j)
            end do 
        end do
    end do
    
    !出力
    open(1, file = 'Bmatrix.dat', status = 'replace')
        do k = 1, ecount
            do l = 1, ipn
                write(1, '(A, I2, A, I2)') "element_", k, "-", l
                do i = 1, 6
                    if(mode == 1) then
                        write(1, '(24F15.7)') (Bmatrix(k, l, i, j), j = 1, scount*dim)
                    else if(mode == 2) then
                        write(1, '(33F15.7)') (Bmatrix(k, l, i, j), j = 1, scount*dim)
                    end if
                end do
            end do
        end do
    close(1)
    
    print *, 'made Bmatrix(output file : "Bmatrix.dat")'

end subroutine