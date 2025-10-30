subroutine make_Beta

    use set_parameter
    implicit none
    
    integer n,i,j,k,l
    double precision Beta_sum     !Ep_ac計算用
    
    !配列初期化
    hf= 0.0d0
     S_dev = 0.0d0
    S_equ = 0.0d0
     Beta_sum = 0.0d0
  if(step == 1) then
        Beta = 0.0d0        !粘塑性関数（解析開始前は加工硬化が起こらないので0）
        Ep_acc = 0.0d0      !蓄積された粘塑性ひずみ（accumulated visco-plastic strain）
 
  
    else
        !蓄積された粘塑性ひずみ
        !εp(ti) = εp(ti-1) + (2/3 * (βij(ti-1)*βij(ti-1)))^0.5 * dt
 do i=1,scount
   do k=1,ecount
   Ep_acc(k,i) = Ep_acc(k,i) + ((Beta(k,i,1)*Beta(k,i,1) + Beta(k,i,2)*Beta(k,i,2) + Beta(k,i,3)*Beta(k,i,3) &
                              + 0.50d0 * (Beta(k,i,4)*Beta(k,i,4) + Beta(k,i,5)*Beta(k,i,5) &
                              + Beta(k,i,6)*Beta(k,i,6)) * 2.0d0/3.0d0) ** 0.50d0) * dt
   end do
  end do            
  !要素の粘塑性ひずみ
   do k=1,ecount
   Ep_el(k) = sum(Ep_acc(k,:))/8.0d0
   end do
  
        !硬化関数 g(εp) = A * (εp)**B + C
         do n = 1, ecount
 if(material_number(n) == 1) then
            do i=1,scount
                hf(n,i) = hf_p1 * (Ep_acc(n,i) ** hf_p2) + hf_p3 
        end do
        
        else if(material_number(n) == 2) then
   
            do i=1,scount
                hf(n,i) = hf2_p1 * (Ep_acc(n,i) ** hf2_p2) + hf2_p3
               
        end do
        end if
        end do
        

    
        !偏差応力 sij(ti) = σ(ti-1) - 1/3 * (δij*σkk(ti-1))
        !汎用性なし、要検討
        do n=1,ecount
            do i=1,scount
                 S_dev(n,i,1) =  microS_gauss(n,i,1) - ( microS_gauss(n,i,1) +  microS_gauss(n,i,2) +  microS_gauss(n,i,3)) / 3.0d0
                 S_dev(n,i,2) =  microS_gauss(n,i,2) - ( microS_gauss(n,i,1) +  microS_gauss(n,i,2) +  microS_gauss(n,i,3)) / 3.0d0
                 S_dev(n,i,3) =  microS_gauss(n,i,3) - ( microS_gauss(n,i,1) +  microS_gauss(n,i,2) +  microS_gauss(n,i,3)) / 3.0d0
                 S_dev(n,i,4) =  microS_gauss(n,i,4)
                 S_dev(n,i,5) =  microS_gauss(n,i,5)
                 S_dev(n,i,6) =  microS_gauss(n,i,6)
            end do
        end do
        
        !相当応力σeq(ti) = (3/2 * sij(ti) * sij(ti))**1/2
        !偏差応力sのせん断成分は2倍で計算
        do n=1,ecount
            do i=1,scount
                do j=1,dir/2
                    S_equ(n,i) = S_equ(n,i) + ( S_dev(n,i,j)**2.0d0)
                end do
                do j=dir/2+1,dir
                    S_equ(n,i) = S_equ(n,i) + (( S_dev(n,i,j)*2.0d0)**2)
                end do
            end do
        end do
    
        do n=1,ecount
            do i=1,scount
                S_equ(n,i) = (1.50d0 * S_equ(n,i)) ** 0.50d0
            end do
        end do
        
        !粘塑性関数β(ti) = 3/2 * ε0 * (σe(ti) / g(εp(ti)) )^n * sij(ti) / σeq(ti)
        !参照ひずみ速度ε0は与える
        do k = 1, ecount
            do l = 1, ipn
                do i = 1, 3
                    Beta(k, l, i) = 1.50d0 * e_rate(material_number(k)) * (S_equ(k, l) / hf(k, l)) ** nonl * S_dev(k, l, i) / S_equ(k, l)
                end do
                
                do i = 4, 6     !テンソルせん断ひずみから工学せん断ひずみに戻す（予想）
                    Beta(k, l, i) = 3.0d0 * e_rate(material_number(k)) * (S_equ(k, l) / hf(k, l)) ** nonl * S_dev(k, l, i) / S_equ(k, l)
                end do
            end do
        end do
        
    end if
    
    
    
    !粘塑性関数の出力
    open(10,file='Beta.dat',status='replace')
        write(FMT, '("("I0"f20.10)")') dir
        do n=1,ecount
            write(10,'(A,I5)')"beta_",n
            do i=1,scount
                write(10,FMT) (beta(n,i,j),j=1,dir)
            end do
            write(10,'(/s)')
        end do
    close(10)
    open(10,file='Beta_ele.dat',status='replace')
        write(FMT, '("("I0"f20.10)")') dir
        do n=1,ecount
            write(10,'(A,I20)')"beta_",n
          
                write(10,FMT) (beta(n,1,j),j=1,dir)
           
            write(10,'(/s)')
        end do
    close(10)
    !硬化関数の出力
    open(11,file='hf.dat',status='replace')
        write(FMT, '("("I0"f20.10)")') scount
        do n=1,ecount
            write(11,'(A,I5)')"hf_",n
            do i=1,scount
                write(11,FMT) (hf(n,i),j=1,scount)
            end do
            write(11,'(/s)')
        end do
    close(11)
    
    !偏差応力の出力
    open(12,file=' S_dev.dat',status='replace')
        write(FMT, '("("I0"f20.10)")') dir
        do n=1,ecount
            write(12,'(A,I5)')"beta_",n
            do i=1,scount
                write(12,FMT) ( S_dev(n,i,j),j=1,dir)
            end do
            write(12,'(/s)')
        end do
    close(12)
    
    !相当応力の出力
    open(13,file='Stress_e.dat',status='replace')
        write(FMT, '("("I0"f20.10)")') scount
        do n=1,ecount
            write(13,'(A,I5)')"beta_",n
            write(13,FMT) ( S_equ(n,i),i=1,scount)
            write(13,'(/s)')
        end do
    close(13)
    !粘塑性ひずみの出力
        open(13,file='Ep_el.dat',status='replace')
        write(FMT, '("("I0"f20.10)")') scount
        do n=1,ecount
            write(13,FMT) (n,Ep_el(n))
        end do
    close(13)
!   do n=1,ecount
!        Ep_el(n)=(Ep_acc(n,1)+Ep_acc(n,2)+Ep_acc(n,3)+Ep_acc(n,4)+Ep_acc(n,5)+Ep_acc(n,6)+Ep_acc(n,7)+Ep_acc(n,8))/8
!        end do
end subroutine
