subroutine make_Kmatrix

    use set_parameter
    
    implicit none
    
    integer i, j, k, i1, i2, j1, j2
    
    !初期化
    Kmatrix = 0.0d0
    Kmat_LU = 0.0d0
    
    !eKmatrixの作成
        call make_Hmatrix
        call make_invJmatrix
        call make_cartesian
        call make_Bmatrix
        call make_eKmatrix

     ![eK]を足し合わせることでKmatrixを作成
     do k = 1, ecount
        do i1 = 1, ipn
            do i2 = 1, ipn
                do j1 = 1, dim
                    do j2 = 1, dim
                        Kmatrix((element(k, i1)-1)*3+j1, (element(k, i2)-1)*3+j2) = Kmatrix((element(k, i1)-1)*3+j1, (element(k, i2)-1)*3+j2) + eKmatrix(k, (i1-1)*3+j1, (i2-1)*3+j2)
                    end do
                end do
            end do
        end do
     end do
     
     !出力
     open(1, file = 'Kmatrix.dat', status = 'replace')
        do i = 1, ncount*dim
            write(1, '(192F15.7)') (Kmatrix(i, j), j = 1, ncount*dim)
        end do
     close(1)
     
     print *, 'made Kmatrix (output file : "Kmatrix.dat")'
          
end subroutine