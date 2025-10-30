subroutine make_Dmatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, e
    
    !Dmatrixの初期化
    Dmatrix = 0.0d0
    
    !等方性材料の材料構成則作成
    do k = 1, ecount
        do i = 1, 3
            do j = 1, 3
                if(i == j) then
                    Dmatrix(k, i, j) = young(k) * (1.0d0 - poisson(k)) / (1.0d0 + poisson(k)) / (1.0d0 -2.0d0 * poisson(k))
                else
                    Dmatrix(k, i, j) = young(k) * poisson(k) / (1.0d0 + poisson(k)) / (1.0d0 - 2.0d0 * poisson(k))
                end if
            end do
        end do
        
        do i = 4, 6
            Dmatrix(k, i, i) = young(k) / (1.0d0 + poisson(k)) / 2.0d0
        end do
    end do
    
    !異方性材料（繊維)の材料構成則作成
    call make_Smatrix
    
    !繊維要素のDマトリックスに繊維の構成則(invSmatrix)を代入
    do k = 1, fiber_unit
        do i = 1, dir
            do j = 1, dir
                Dmatrix(fiber_number(k), i, j) = invSmatrix(material_number(fiber_number(k))-2, i, j)
            end do
        end do
    end do
    
    open(1, file = 'Dmatrix.dat', status = 'replace')
        do k = 1, ecount
            write(1, '("element_", I3)') k
            do i = 1, dir
                write(1, '(6F15.7)') (Dmatrix(k, i, j), j = 1, dir)
            end do
        end do
    close(1)
    
    print *, 'made Dmatrix'

end subroutine