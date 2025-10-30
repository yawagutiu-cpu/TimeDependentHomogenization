subroutine solve_micro

    use set_parameter
    
    implicit none
    
    integer i, j, k, l, m
    
    !初期化
    macroE_tensor = 0.0d0       !巨視的ひずみテンソル
    microS_element = 0.0d0      !微視的応力(要素ごと)
    microS_node = 0.0d0         !微視的応力(節点ごと)
    microS_vm_e = 0.0d0         !微視的なフォンミーゼス応力(要素ごと)
    microS_vm_g = 0.0d0         !微視的なフォンミーゼス応力(積分点ごと)
    microS_vm_n = 0.0d0         !微視的なフォンミーゼス応力(節点ごと)
    n_counter = 0.0d0           !重複している節点の数を数える
    microS_velocity = 0.0d0     !微視的応力速度
    u_sharp_velocity = 0.0d0    !擾乱変位速度
    node_after = 0.0d0          !変形後の座標
     microU = 0.0d0

    !巨視的ひずみベクトルをテンソルにする
    !工学せん断ひずみをテンソルせん断ひずみに変換
    macroE_tensor(1, 1) = macroE(1)
    macroE_tensor(2, 2) = macroE(2)
    macroE_tensor(3, 3) = macroE(3)
    macroE_tensor(1, 2) = macroE(4) * 5.0d-1
    macroE_tensor(2, 3) = macroE(5) * 5.0d-1
    macroE_tensor(3, 1) = macroE(6) * 5.0d-1
    macroE_tensor(2, 1) = macroE_tensor(1, 2)
    macroE_tensor(3, 2) = macroE_tensor(2, 3)
    macroE_tensor(1, 3) = macroE_tensor(3, 1)
    
    !微視的応力速度の計算
    !microS_velocity = [eA] * macroE_velocity - {eR}
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dir
                do m = 1, dir
                    microS_velocity(k, l, i) = microS_velocity(k, l, i) + eAmatrix(k, l, i, m) * macroE_velocity(m)
                end do
        
                microS_velocity(k, l, i) = microS_velocity(k, l, i) - eRvector(k, l, i)
            end do
        end do
    end do
    
    !微視的応力の計算
    !microS(ti) = microS(t(i-1)) + microS_velocity(ti) * dt
    !積分点ごと，要素ごと

    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dir
                microS_gauss(k, l, i) = microS_gauss(k, l, i) + microS_velocity(k, l, i) * dt   
                microS_element(k, i) = microS_element(k, i) + microS_gauss(k, l, i) / ipn  
              
           
            end do
        end do
    end do
        
    !節点ごと
    !重複している節点の数を数える
    do k = 1, ecount
        do l = 1, ipn
            n_counter(element(k, l)) = n_counter(element(k, l)) + 1 
        end do
    end do
    
    !重複を許して足し合わせる
    do k = 1, ecount
        do l = 1, ipn
            do i = 1, dir
                microS_node(element(k, l), i) = microS_node(element(k, l), i) + microS_gauss(k, l, i)
            end do
        end do
    end do
    
    !重複している数で割ることで平均的な節点ごとの微視的応力を得る((ecount * scount = 216)→(ncount = 64))
    do i = 1, ncount
        do j = 1, dir
            microS_node(i, j) = microS_node(i, j) / n_counter(i)
        end do
    end do
    
    !求めた各微視的応力についてフォンミーゼス応力を求める
    do k = 1, ecount
        do l = 1, ipn
            microS_vm_e(k) = sqrt(5.0d-1 * ((microS_element(k, 1) - microS_element(k, 2)) ** 2.0d0 + (microS_element(k, 2) - microS_element(k, 3)) ** 2.0d0 + (microS_element(k, 3) - microS_element(k, 1)) ** 2.0d0 + 6.0d0 * (microS_element(k, 4) ** 2.0d0 + microS_element(k, 5) ** 2.0d0 + microS_element(k, 6) ** 2.0d0)))
            microS_vm_g(k, l) = sqrt(5.0d-1 * ((microS_gauss(k, l, 1) - microS_gauss(k, l, 2)) ** 2.0d0 + (microS_gauss(k, l, 2) - microS_gauss(k, l, 3)) ** 2.0d0 + (microS_gauss(k, l, 3) - microS_gauss(k, l, 1)) ** 2.0d0 + 6.0d0 * (microS_gauss(k, l, 4) ** 2.0d0 + microS_gauss(k, l, 5) ** 2.0d0 + microS_gauss(k, l, 6) ** 2.0d0)))
        end do
    end do
    
    do i = 1, ncount
        microS_vm_n(i) = sqrt(5.0d-1 * ((microS_node(i, 1) - microS_node(i, 2)) ** 2.0d0 + (microS_node(i, 2) - microS_node(i, 3)) ** 2.0d0 + (microS_node(i, 3) - microS_node(i, 1)) ** 2.0d0 + 6.0d0 * (microS_node(i, 4) ** 2.0d0 + microS_node(i, 5) ** 2.0d0 + microS_node(i, 6) ** 2.0d0)))
    end do
    
    !擾乱変位速度
    !u#_velocity = [χ] * macroE_velocity + {P}
    do i = 1, ncount
        do j = 1, dim
            do m = 1, dir
                u_sharp_velocity(i, j) = u_sharp_velocity(i, j) + Xmatrix((i-1)*dim + j, m) * macroE_velocity(m)
            end do
            
            u_sharp_velocity(i, j) = u_sharp_velocity(i, j) + Pvector((i-1)*dim + j)
        end do
    end do
    
    !擾乱変位の計算
    !u#(ti) = u#(t(i-1)) + u#_velocity(ti) * dt
    do i = 1, ncount
        do j = 1, dim
            u_sharp(i, j) = u_sharp(i, j) + u_sharp_velocity(i, j) * dt
        end do
    end do
    
    !微視的変位の計算
    !先生博論式(2.4)
    do i = 1, ncount
        do j = 1, dim
            do k = 1, dim
                microU(i, j) = microU(i, j) + macroE_tensor(j, k) * node(i, k)
            end do
            microU(i, j) = microU(i, j) + u_sharp(i, j)
        end do
    end do
    
    !変形後の節点座標
    do i = 1, ncount
        do j = 1, dim
            node_after(i, j) = node(i, j) + microU(i, j) * eps
        end do
    end do
    
end subroutine