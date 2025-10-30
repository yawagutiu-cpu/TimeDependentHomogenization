subroutine input

    use set_parameter
    
    implicit none
    
    integer i, j, k, e
    
    !基本情報の読み込み
    open(1, file = 'basic_info.txt', form = 'formatted')
        read(1, *) LABEL, ecount            !要素数
        read(1, *) LABEL, ncount            !節点数
        read(1, *) LABEL, dim               !次元数
        read(1, *) LABEL, eps               !擾乱の倍数
        read(1, *) LABEL, lambda            !ペナルティ数
        read(1, *) LABEL, pair              !ペアリング数
        read(1, *) LABEL, mode              !CE(1) or NCE(2)
        read(1, *) LABEL, Vol

        close(1)
          open(2, file = 'basic_info.dat', form = 'formatted')
            write(2, *) LABEL, ecount            !要素数
         write(2, *) LABEL, ncount            !節点数
         write(2, *) LABEL, dim               !次元数
   write(2, *) LABEL, eps               !擾乱の倍数
     write(2, *) LABEL, lambda            !ペナルティ数
     write(2, *) LABEL, pair              !ペアリング数
       write(2, *) LABEL, mode              !CE(1) or NCE(2)
         write(2, *) LABEL, Vol
      close(2)
    
    dir = dim + dim * (dim - 1) / 2         !3D: dir=6, 2D: dir=3
    ipn = 4 * (dim - 1)                     !ipn: Isoparametric Number 
                                            !3D: ipn=8, 2D: ipn=4
    if(mode == 1) then
        scount = ipn                        !形状関数の数（適合）
    else if(mode == 2) then
        scount = ipn + dim                  !形状関数の数（非適合）
    end if
    
    print *,ipn
    
    !基本情報を基に各種マトリックスにメモリを割り当てる
    call allocate_memory
    
    !初期情報(要素，節点)の読み込み
    open(2, file = 'element.txt', form = 'formatted')
        do i = 1, ecount
            read(2, *) k, (element(i, j), j = 1, ipn)
        end do
    close(2)
    
    open(3, file = 'node.txt', form = 'formatted')
        do i = 1, ncount
            read(3, *) k, (node(i, j), j = 1, dim)
        end do
    close(3)
    
    !ユニットセルの体積計算
    !x_temp = 0.0d0
    !y_temp = 0.0d0
    !z_temp = 0.0d0
    
    !do k = 1, ncount
     !   if(x_temp < node(k, 1)) then
      !      x_temp = node(k, 1)
       ! end if
        !if(y_temp < node(k, 2)) then
         !   y_temp = node(k, 2)
        !end if
        !if(z_temp < node(k, 3)) then
         !   z_temp = node(k, 3)
        !end if
    !end do
    
    !Vol = x_temp * y_temp * z_temp
  
    
    !要素材質決定配列(1:アルミ合金，2:エポキシ，3:1軸方向炭素繊維，4:２軸方向炭素繊維，5:３軸方向炭素繊維)
    open(4, file = 'element_type.txt', form = 'formatted')
        do i = 1, ecount
            read(4, *) k, material_number(i)
        end do
    close(4)
    
    !可変材料定数(アルミ合金)の読み込み
    open(5, file = 'material_property.txt', form = 'formatted')
    open(6,file='material_info.dat',form='formatted')
        read(5, *) LABEL, young1, poisson1
        read(5, *) LABEL, young2, poisson2
        read(5, *) LABEL, Ef_1, Ef_3, poissonf_31, poissonf_12, Gf_31, nonl, hf_p1, hf_p2, hf_p3
        read(5, *) LABEL, Ef2_1, Ef2_3, poissonf2_31, poissonf2_12, Gf2_31, nonl2, hf2_p1, hf2_p2, hf2_p3
                          !(1.55d4, 2.40d5, 0.28d0, 0.49d0, 2.47d4)
                           write(6,*) '1のヤング率',young1,poisson1
                           write(6,*) '2のヤング率',young2,poisson2
                            write(6,*) '累乗数',nonl
                            write(6,*)  '1の硬化関数',hf_p1, hf_p2, hf_p3
                            write(6,*) '2の硬化関数', hf2_p1, hf2_p2, hf2_p3
     close(6)
    close(5)
    
    !等方性材料(軟鋼 or エポキシ)要素の材料定数の決定
    do i = 1, ecount
        if(material_number(i) == 1) then
            young(i) = young1
            poisson(i) = poisson1
        else if(material_number(i) == 2) then
            young(i) = young2
            poisson(i) = poisson2
        else
            young(i) = 0.0d0    !炭素繊維については別で計算する
            poisson(i) = 0.0d0
        end if
    end do
    
    !異方性材料(炭素繊維)要素情報の整理
    do i = 1, ecount
        if(material_number(i) == 3 .or. material_number(i) == 4 .or. material_number(i) == 5) then
            fiber_unit = fiber_unit + 1 !炭素繊維要素の数を数える
        end if
    end do
    
    allocate (fiber_number(fiber_unit))
    
    j = 1
    do i = 1, ecount
        if(material_number(i) == 3 .or. material_number(i) == 4 .or. material_number(i) == 5) then
            fiber_number(j) = i         !炭素繊維要素の要素番号を保存
            j = j + 1
        end if
    end do
    
    !ひずみ制御に関する設定
    e_rate = 0.0d0
    
    open(6, file = 'strain_control_info.txt', form = 'formatted')
    open(8,file='steain.txt',form='formatted')
        read(6, *) LABEL, dire             !変位速度の方向
        read(6, *) LABEL, dt               !時間の刻み幅
        read(6, *) LABEL, disp_end         !変位目標値
        read(6, *) LABEL, Umat_speed       !変位速度
        read(6, *) LABEL, strain_end       !ひずみ目標値
        read(6, *) LABEL, strain_speed     !与えるひずみ速度
        read(6, *) LABEL, e_rate(1)           !参照ひずみ速度
        write(8,*) dire
         write(8,*) dt
          write(8,*) disp_end
           write(8,*) Umat_speed
            write(8,*) strain_end
             write(8,*) strain_speed
              write(8,*) e_rate
        close(8)
    close(6)
    do i=1,5
   e_rate(i)=e_rate(1)
   end do
       e_rate(2)= e_rate(1)    
    !出力する要素・積分点
    open(7, file = 'check_info.txt', form = 'formatted')
        read(7, *) LABEL, check_e
        read(7, *) LABEL, check_g
    close(7)
    
   ! ペアリングする節点の指定
    open(8, file = 'penalty.txt', form = 'formatted')
        do i = 1, pair
            read(8, *) k, penalty_number(i, 1), penalty_number(i, 2)
        end do
    close(8)    
    
    print *, 'completed input of parameter' 
    
end subroutine