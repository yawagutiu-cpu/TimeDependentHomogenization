subroutine solve_macro

    use set_parameter
    
    implicit none
    
    integer i, m
    
    !初期化
    macroE_velocity = 0.0d0    !巨視的ひずみ速度
    macroS_velocity = 0.0d0    !巨視的応力速度
    macroE_vm = 0.0d0       !フォンミーゼス(Von Mises)ひずみ
    macroS_vm = 0.0d0       !フォンミーゼス(Von Mises)応力
    
    !巨視的ひずみ速度を与える
    !direで指定した方向に一定のひずみ速度strain_speedを与える
    macroE_velocity(dire) = strain_speed
    
    !巨視的応力速度の計算
    !巨視的応力速度のつり合い式：macroS_velocity = <A>*macroE_velocity - <R>
    !単軸引張条件とすると引張方向以外の応力速度成分はゼロになる
    !例えば3軸方向にひずみ速度を与えたとすると, 次のように考えられる
    !--------------------------------------------------------------------------------------------
    !  |  0  |         |    Ev11    |           |    Ev11    |            |  0  |
    !  |  0  |         |    Ev22    |           |    Ev22    |            |  0  |
    !  |Σv33| = <A> * |strain_speed| - <R> ⇔　|strain_speed| = <A>^-1 * |Σv33| + <A>^-1 * <R>
    !  |  0  |         |    Ev12    |           |    Ev12    |            |  0  | 
    !  |  0  |         |    Ev23    |           |    Ev23    |            |  0  |
    !  |  0  |         |    Ev31    |           |    Ev31    |            |  0  |
    !--------------------------------------------------------------------------------------------   
    !33成分だけ取り出せば,
    !strain_speed = a^-1_33 * Σv33 + a^-1_31 * r1 + … + a^-1_36 * r6
    !よって，巨視的応力速度は次のように求まる
    !Σv33 = (strain_speed - a^-1_31 * r1 - … - a^-1_36 * r6) / a^-1_33
    
    macroS_velocity(dire) = macroE_velocity(dire) / invAmatrix(dire, dire)
    
    do m = 1, dir
        macroS_velocity(dire) = macroS_velocity(dire) - invAmatrix(dire, m) * Rvector(m) / invAmatrix(dire, dire)
    end do
    
    !巨視的ひずみ速度の計算
    !macroE_velocity = <A>^-1 * macroS_velocity + <A>^-1 * <R>
 
    do i = 1, dir
        if(i /= dire) then
            do m = 1, dir
                macroE_velocity(i) = macroE_velocity(i) + invAmatrix(i, m) * (macroS_velocity(m) + Rvector(m))
          
                   
            end do
        end if
    end do

    
    !巨視的応力・ひずみの計算
    !macroS,E(ti) = macroS,E(t(i-1)) + macroS,E_velocity(ti) * dt
    do i = 1, dir
        macroS(i) = macroS(i) + macroS_velocity(i) * dt
        macroE(i) = macroE(i) + macroE_velocity(i) * dt
    end do
    
    !フォンミーゼス応力・ひずみの計算
    macroS_vm = sqrt(5.0d-1 * ((macroS(1) - macroS(2)) ** 2.0d0 + (macroS(2) - macroS(3)) ** 2.0d0 + (macroS(3) - macroS(1)) ** 2.0d0 + 6.0d0 * (macroS(4) ** 2.0d0 + macroS(5) ** 2.0d0 + macroS(6) ** 2.0d0)))
    
    macroE_vm = sqrt(5.0d-1 * ((macroE(1) - macroE(2)) ** 2.0d0 + (macroE(2) - macroE(3)) ** 2.0d0 + (macroE(3) - macroE(1)) ** 2.0d0 + 6.0d0 * (macroE(4) ** 2.0d0 + macroE(5) ** 2.0d0 + macroE(6) ** 2.0d0)))

end subroutine