!  main.f90 
!****************************************************************************
!
!  プログラム: Time_Dependent_Homogenization
!
!  目的:  弾-粘塑性を考慮した時間依存均質化理論の構築(ひずみ制御)
!
!  入力： basic_info.txt（要素数・節点数・積分点数・次元数・ユニットセルYの体積・ペナルティ数・ペアリング数）
!         node.txt, element.txt（節点情報, 要素情報）
!         material_number.txt（各要素の材質）
!         material_property.txt（材料定数（軟鋼，エポキシ，炭素繊維））
!         penalty.txt（ペナルティ法においてペアリングする節点情報）
!
!  出力：Dmatrix.dat（応力-ひずみマトリックス）
!        Kmatrix.dat（全体剛性マトリックス）
!        Fmatrix.dat
!        Xmatrix.dat（特性変位マトリックス）
!        Amatrix.dat（巨視的弾性剛性テンソル）
!        Macro_TD.dat（巨視的ひずみ・応力の時間変化）
!        output_***.dat（時間依存ステップ終了後の情報）
!
!****************************************************************************

program main

    use set_parameter

    implicit none

    call input
    call set_coordinates
    call make_Dmatrix
    call make_Kmatrix
    call boundary_condition
    call make_Xmatrix
    call make_Amatrix
    call Time_Dependent
        call output_paraview
    print *, 'program completed'
    
end program 
