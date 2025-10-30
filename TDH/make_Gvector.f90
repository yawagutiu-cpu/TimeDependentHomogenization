subroutine make_Gvector

    use set_parameter
    
    implicit none
    
    integer i, j, k, l, m
    
    !初期化
    eGvector = 0.0d0
    eGvec_c = 0.0d0
    eGvec_n = 0.0d0
    eKcn_inveKnn_eGn = 0.0d0
    Gvector = 0.0d0
    
    !eGvectorの計算
    !{eG} = [B]t[D] * {β(ti)} * det_J
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, scount*dim
                do m = 1, dir
                    eGvector(k, i) = eGvector(k, i) + tBmatDmat(k, l, i, m) * Beta(k, l, m) * det_J(k, l)
                end do
            end do
        end do
    end do
    
    !非適合要素を使用する場合の処理**************************************************************************************************
    if(mode == 2) then
        !関連しているモード（適合か非適合）ごとにベクトルを下記のように分ける
        !------------------------------
        !              |              |
        !              |     eGc      |
        !              |    (24*1)    |
        ! eGvector =   |              |
        !              |--------------|
        !              |     eGn      |
        !              |    (9*1)     |
        !------------------------------    
        do k = 1, ecount
            do i = 1, ipn*dim
                eGvec_c(k, i) = eGvector(k, i)
            end do
            do i = 1, dim*dim
                eGvec_n(k, i) = eGvector(k, ipn*dim + i)
            end do
        end do
        
        !eKmatrixの時と同様に非適合モードに対応する変位場（u9,u10,u11,v9,v10,v11,w9,w10,w11）に関する式を消去する
        ![eKcn]*[eKnn]^-1*{eGn}の計算
        do k = 1, ecount
            do i = 1, ipn*dim
                do m = 1, dim*dim
                    eKcn_inveKnn_eGn(k, i) = eKcn_inveKnn_eGn(k, i) + eKcn_inveKnn(k, i, m) * eGvec_n(k, m)
                end do
            end do
        end do
        
        !便宜上eGvectorに上書きする
        eGvector = 0.0d0
        
        do k = 1, ecount
            do i = 1, ipn*dim
                eGvector(k, i) = eGvec_c(k, i) - eKcn_inveKnn_eGn(k, i)
            end do
        end do
    end if
    !********************************************************************************************************************************
    
    !Gvectorの作成
    do k = 1, ecount
        do i = 1, ipn
            do j = 1, dim
                Gvector((element(k, i)-1)*dim + j) = Gvector((element(k, i)-1)*dim + j) + eGvector(k, (i-1)*dim + j)
            end do
        end do
    end do

end subroutine