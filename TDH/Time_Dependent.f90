subroutine Time_Dependent

    use set_parameter
    
    implicit none
    
    integer i, j, k, l
    

    
    !時間依存ステップの初期設定
    swich_end = 0               !最終ステップスイッチ
    step = 1                    !ステップ数
    ti = 0.0d0                  !時間
    strain_giv = strain_speed*dt  !与えるひずみの初期値
    macroE = 0.0d0
    macroS = 0.0d0
    microS_gauss = 0.0d0
    u_sharp = 0.0d0
    microU = 0.0d0
    
     !各ステップで求めた値をその都度記録する
    open(1, file = 'Macro_TD.dat', status = 'replace')
    open(2, file = 'VonMises_Macro_TD.dat', status = 'replace')
        write(1, FMT = '(a)', advance = 'no') '         ti          macroE(1)      macroS(1)      macroE(2)      macroS(2)      macroE(3)      macroS(3)      macroE(4)      macroS(4)      macroE(5)      macroS(5)      macroE(6)      macroS(6)' 
        write(1, *)
      
        write(2, FMT = '(a)', advance = 'no') '         ti          macroE_vm      macroS_vm'
        write(2, *)
 
        !与えるひずみが設定値になるまで計算を繰り返す
        do while(strain_giv < strain_end)
            
            !前のループの計算結果を出力
            write(1, '(13F15.7)') ti, ((macroE(i), macroS(i)), i = 1, dir)
         
            write(2, '(3F15.7)') ti, macroE_vm, macroS_vm
  
            call make_Beta
            call make_Gvector
            call make_Pvector
            call make_Rvector
            call solve_macro
            call solve_micro
           
            step = step + 1
            ti = ti + dt
            strain_giv = strain_giv + strain_speed * dt
        end do
        
        write(1, '(13F15.7)') ti, ((macroE(i), macroS(i)), i = 1, dir)
             
        write(2, '(3F15.7)') ti, macroE_vm, macroS_vm
                    
        !最終ステップの計算***************************************************************************************************************************************************************************************************************************
        swich_end = 1
        
        print *, "step_end : ", step
        print *, "ti_end : ", ti+dt
        print *, "strain_giv_end : ", strain_giv
        
        call make_Beta
        call make_Gvector
        call make_Pvector
        call make_Rvector
        call solve_macro
        call solve_micro
        
        write(1, '(13F15.7)') ti+dt, ((macroE(i), macroS(i)), i = 1, dir)
        
        write(2, '(3F15.7)') ti+dt, macroE_vm, macroS_vm
              
    close(1)
    close(2)
      close(3)
    
    !最終計算結果の出力*******************************************************************************************************************************************************************************************************************************
            open(3, file = 'Ep_acc.dat', status = 'replace')
        do k = 1, ecount
    do l = 1, ipn
        write(3, *) k, l, Ep_acc(k, l)
    end do
end do
close(3)
      open(11,file='el_hizumi.dat',status='replace')
        do k=1,ecount
            write(11,'(I5,E15.2e3)')k,Ep_el(k)
            end do
            write(11,'(/s)')
close(11)
      !粘塑性関数，硬化関数，偏差応力，相当応力
  !  open(3, file = 'output_Beta.dat', status = 'replace')
   ! open(4, file = 'output_Hardening_Function.dat', status = 'replace')
    !open(5, file = 'output_Deviatoric_Stress.dat', status = 'replace')
    !open(6, file = 'output_Equivalent_Stress.dat', status = 'replace')
     !   open(7, file = ' e_rate(material_number).dat', status = 'replace')
      !  do k = 1, ecount
       !     write(3, '("element_", I4)') k
        !    write(3, FMT = '(a)', advance = 'no') '            ipn_1               ipn_2               ipn_3               ipn_4               ipn_5               ipn_6               ipn_7               ipn_8'
         !   write(3, *)
          !  write(4, '("element_", I4)') k
           ! write(4, FMT = '(a)', advance = 'no') '            ipn_1               ipn_2               ipn_3               ipn_4               ipn_5               ipn_6               ipn_7               ipn_8'
            !write(4, *)
            !write(5, '("element_", I4)') k
            !write(5, FMT = '(a)', advance = 'no') '            ipn_1               ipn_2               ipn_3               ipn_4               ipn_5               ipn_6               ipn_7               ipn_8'
            !write(5, *)
            !write(6, '("element_", I4)') k
            !write(6, FMT = '(a)', advance = 'no') '            ipn_1               ipn_2               ipn_3               ipn_4               ipn_5               ipn_6               ipn_7               ipn_8'
            !write(6, *)
            !write(7,*) e_rate(material_number(k))
            !do i = 1, dir
             !   write(3, '(8F20.10)') (Beta(k, l, i), l = 1, ipn)
              !  write(5, '(8F20.10)') (S_dev(k, l, i), l = 1, ipn) 
           ! end do
            !write(4, '(8F20.10)') (hf(k, l), l = 1, ipn)
            !write(6, '(8F20.10)') (S_equ(k, l), l = 1, ipn)        
        !end do
    !close(3)
    !close(4)
    !close(5)
    !close(6)
    !close(7)
    
    !ePvector, Pvector
    open(9, file = 'output_ePvector.dat', status = 'replace')
    open(10, file = 'output_Pvector.dat', status = 'replace')
        do k = 1, ecount
            write(9, '("element_", I4)') k
            do i = 1, scount*dim
                write(9, '(F20.10)') ePvector(k, i)
            end do
        end do
        do i = 1, ncount*dim
            write(10, '(I3, F20.10)') i, Pvector(i)
        end do
    close(9)
    close(10)
    
    !eRvector, Rvector
    open(11, file = 'output_eRvector.dat', status = 'replace')
    open(12, file = 'output_Rvector.dat', status = 'replace')
        do k = 1, ecount
            write(11, '("element_", I4)') k
            write(11, FMT = '(a)', advance = 'no') '            ipn_1               ipn_2               ipn_3               ipn_4               ipn_5               ipn_6               ipn_7               ipn_8'
            write(11, *)
                do i = 1, dir
                    write(11, '(8F20.10)') (eRvector(k, l, i), l = 1, ipn)
                end do
        end do
        do i = 1, dir
            write(12, '(F20.10)') Rvector(i)
        end do
    close(11)
    close(12)
    
    !巨視的ひずみテンソル
    open(13, file = 'output_macroE_tensor.dat', status = 'replace')
        do i = 1, dim
            write(13, '(3F20.10)') (macroE_tensor(i, j), j = 1, dim)
        end do
    close(13)
    
    !微視的応力
    open(14, file = 'output_microS_gauss.dat', status = 'replace')
    open(15, file = 'output_microS_element.dat', status = 'replace')
    open(16, file = 'output_microS_node.dat', status = 'replace')
        do k = 1, ecount
            write(14, '("element_", I4)') k
            write(14, FMT = '(a)', advance = 'no') '            ipn_1               ipn_2               ipn_3               ipn_4               ipn_5               ipn_6               ipn_7               ipn_8'
            write(14, *)
            write(15, '("element_", I4)') k
            do i = 1, dir
                write(14, '(8F20.10)') (microS_gauss(k, l, i), l = 1, ipn)
                write(15, '(F20.10)') microS_element(k, i)
            end do
        end do
        do i = 1, ncount
            write(16, '("node_", I4)') i
            do j = 1, dir
                write(16, '(F20.10)') microS_node(i, j)
            end do
        end do
    close(14)
    close(15)
    close(16)
    open(15, file = 'output_microS_ele-z.dat', status = 'replace')
     do k = 1, ecount
      write(15, '("element_", I4)') k
        write(15, '(F20.10)') microS_element(k,3)
        end do
        close(15)
            open(15, file = 'output_microS_ele-x.dat', status = 'replace')
     do k = 1, ecount
      write(15, '("element_", I4)') k
        write(15, '(F20.10)') microS_element(k,1)
        end do
        close(15)
            open(15, file = 'output_microS_ele-y.dat', status = 'replace')
     do k = 1, ecount
      write(15, '("element_", I4)') k
        write(15, '(F20.10)') microS_element(k,2)
        end do
        close(15)
            !擾乱変位, 微視的変位
    open(17, file = 'output_u_sharp.dat', status = 'replace')
    open(18, file = 'output_microU.dat', status = 'replace')
        do i = 1, ncount
            write(17, '(I4, 3F20.10)') i, (u_sharp(i, j), j = 1, dim)
            write(18, '(I4, 3F20.10)') i, (microU(i, j), j = 1, dim)
        end do
    close(17)
    close(18)
    
    !変形後の座標
    open(19, file = 'output_node_after.dat', status = 'replace')
        do i = 1, ncount
            write(19, '(I4, 3F20.10)') i, (node_after(i, j), j = 1, dim)
        end do
    close(19)
    
    !微視的なフォンミーゼス応力
    open(20, file = 'output_VonMises_microS_gauss.dat', status = 'replace')
    open(21, file = 'output_VonMises_microS_element.dat', status = 'replace')
    open(22, file = 'output_VonMises_microS_node.dat', status = 'replace')
        do k = 1, ecount
            write(20, '("element_", I4)') k
            write(20, FMT = '(a)', advance = 'no') '            ipn_1               ipn_2               ipn_3               ipn_4               ipn_5               ipn_6               ipn_7               ipn_8'
            write(20, *)
            write(20, '(8F20.10)') (microS_vm_g(k, l), l = 1, ipn)        
        end do
        do i = 1, ecount
            write(21, '(I4, F20.10)') i, microS_vm_e(i)
        end do
        do i = 1, ncount
            write(22, '(I4, F20.10)') i, microS_vm_n(i)
        end do
    close(20)
    close(21)
    close(22)
    
    print *, 'finished Time Dependent Homogenization'

end subroutine