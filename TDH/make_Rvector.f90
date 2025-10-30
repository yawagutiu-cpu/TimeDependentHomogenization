subroutine make_Rvector

    use set_parameter
    
    implicit none
    
    integer i, j, k, l, m, e
    
    !初期化
    DmatBeta = 0.0d0
    DmatBmat = 0.0d0
    eRvector = 0.0d0
    Rvector = 0.0d0
    
    !eRvectorの計算
    !{eR} = [D]*{β} - [D][B]*{P(ti)} 
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dir
                do m = 1, dir
                    DmatBeta(k, l, i) = DmatBeta(k, l, i) + Dmatrix(k, i, m) * Beta(k, l, m)
                end do
            end do
        end do
    end do
    
    do k = 1, ecount
        do l = 1, ipn 
            do i = 1, dir
                do j = 1, scount*dim
                    do m = 1, dir
                        DmatBmat(k, l, i, j) = DmatBmat(k, l, i, j) + Dmatrix(k, i, m) * Bmatrix(k, l, m, j)
                    end do
                end do
            end do
        end do
    end do
    
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dir
                eRvector(k, l, i) = DmatBeta(k, l, i)

                do m = 1, scount*dim
                    eRvector(k, l, i) = eRvector(k, l, i) - DmatBmat(k, l, i, m) * ePvector(k, m)
                end do
            end do
        end do
    end do
       
    !Rvectorの作成
    do i = 1, dir
        do k = 1, ecount
            do l = 1, ipn
                Rvector(i) = Rvector(i) + eRvector(k, l, i) * det_J(k, l)
            end do
        end do
    end do
    
    !均質化
    Rvector = Rvector / Vol

end subroutine