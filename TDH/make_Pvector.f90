subroutine make_Pvector

    use set_parameter
    
    implicit none
    
    integer i, k, l, m, lda, n, ldb, nrhs
    
    !初期化
    Pvector = 0.0d0
    ePvec_c = 0.0d0
    eKnc_ePc = 0.0d0
    ePvec_n = 0.0d0
    ePvector = 0.0d0
    
    !粘塑性に関する特性関数φ(Pvector)の作成
    ![K]{φ}={G}を解く
    !毎計算ステップで求め直さなければならないので事前にLU分解しておいたKmat_LUを利用する
    Pvector = Gvector
    
    lda = ubound(Kmatrix, 1)
    n = ubound(Kmatrix, 2)
    ldb = ubound(Pvector, 1)
    nrhs = 1
    
    call dgetrs(trans, n, nrhs, Kmat_LU, lda, ipiv_k, Pvector, ldb, info)
        
    !ePvectorの作成
    !適合要素を使用する場合の処理*******************************************************************************************************
    if(mode == 1) then
        !PvectorをePvectorに入れ直す
        do k = 1, ecount
            do l = 1, ipn
                do i = 1, dim
                    ePvector(k, (l-1)*dim + i) = Pvector((element(k, l)-1)*dim + i)
                end do
            end do
        end do
    !非適合要素を使用する場合の処理**************************************************************************************************
    else if(mode == 2) then
        !節点ごとになっているPvectorを要素ごとに振り分けることでePvec_cを得る
        do k = 1, ecount
            do l = 1, ipn
                do i = 1, dim
                    ePvec_c(k, (l-1)*dim + i) = Pvector((element(k, l)-1)*dim + i)
                end do
            end do
        end do
        
        !ePvec_cからePvec_nを求める
        !{ePn}=[eKnn]^-1*({eGn}-[eKnc]{ePc})
        do k = 1, ecount
            do i = 1, dim*dim
                do m = 1, ipn*dim
                    eKnc_ePc(k, i) = eKnc_ePc(k, i) + eKmat_nc(k, i, m) * ePvec_c(k, m)
                end do
            end do
        end do
        
        do k = 1, ecount
            do i = 1, dim*dim
                do m = 1, dim*dim
                    ePvec_n(k, i) = ePvec_n(k, i) + inveKmat_nn(k, i, m) * (eGvec_n(k, m) - eKnc_ePc(k, m))
                end do
            end do
        end do
        
        !求めたePvec_c,nを基にePvectorを作成
        !------------------------------
        !              |              |
        !              |     ePc      |
        !              |    (24*1)    |
        ! ePvector =   |              |
        !              |--------------|
        !              |     ePn      |
        !              |    (9*1)     |
        !------------------------------
        do k = 1, ecount
            do i = 1, ipn*dim
                ePvector(k, i) = ePvec_c(k, i)
            end do
            do i = 1, dim*dim
                ePvector(k, ipn*dim + i) = ePvec_n(k, i)
            end do
        end do    
    end if
    !********************************************************************************************************************************

end subroutine