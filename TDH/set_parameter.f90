module set_parameter

    implicit none
    
    !基本情報
    integer ecount     !要素数
    integer ncount     !節点数
    integer scount     !形状関数の数
    integer dim        !次元数
    integer ipn        !1要素の節点数
    integer dir        !dim + dim * (dim - 1) / 2
    
    integer mode       !適合要素か非適合要素か
    
    integer check_e    !出力する要素
    integer check_g    !出力する積分点
        
    double precision x_temp
    double precision y_temp
    double precision z_temp
    double precision Vol  !ユニットセルYの体積
    double precision eps  !擾乱の倍数
    double precision, parameter :: pi = 3.141592654d0        !円周率
    double precision, parameter :: ga = 0.577350269189626d0  !ガウス点
    
    character(len=100) FMT
    character(len=100) LABEL
    
    integer, allocatable, dimension(:,:) :: element
    double precision, allocatable, dimension(:,:) :: node
    double precision, allocatable, dimension(:,:) :: gauss
    double precision, allocatable, dimension(:,:) :: polar
    double precision, allocatable, dimension(:,:) :: det_J

    !ひずみ制御に関する設定値    
    double precision dire             !変位速度の方向
    double precision dt               !時間の刻み幅
    double precision disp_end         !変位目標値
    double precision Umat_speed       !変位速度
    double precision strain_end       !ひずみ目標値
    double precision strain_speed     !与えるひずみ速度
    double precision e_rate(1:5)           !参照ひずみ速度
        
    !境界条件
    double precision lambda   !ペナルティ数(節点を固定するために代入する値)
    integer pair               !ペアリング数(pair=dim*(dim+1)*(dim+1))
    integer, allocatable, dimension(:,:) :: penalty_number
    double precision, allocatable, dimension(:,:) :: x
    
    !材料情報
    double precision, allocatable, dimension(:) :: material_number     !材料番号
    double precision, allocatable, dimension(:) :: young               !要素ごとのヤング率
    double precision, allocatable, dimension(:) :: poisson             !要素ごとのポアソン比
    integer fiber_unit     !炭素繊維要素の数
    integer, allocatable, dimension(:) :: fiber_number  !炭素繊維の要素番号
    
    !等方性材料の材料定数
    double precision young1
    double precision poisson1
    double precision young2
    double precision poisson2
    
    !異方性材料の材料定数(右手座標系x3軸が繊維方向)
    double precision :: Ef_3             !炭素繊維(carbon fiber)繊維方向ヤング率
    double precision :: Ef_1             !炭素繊維半径方向ヤング率(Ef_2=Ef_1)
    double precision :: poissonf_31      !炭素繊維ポアソン比(νf_32=νf_31)
    double precision :: poissonf_12      !炭素繊維ポアソン比
    double precision :: Gf_31            !炭素繊維剛性率(Gf_23=Gf_31)
    double precision :: nonl
    double precision :: hf_p1
    double precision :: hf_p2
    double precision :: hf_p3
    double precision :: Ef2_3             !炭素繊維(carbon fiber)繊維方向ヤング率
    double precision :: Ef2_1             !炭素繊維半径方向ヤング率(Ef_2=Ef_1)
    double precision :: poissonf2_31      !炭素繊維ポアソン比(νf_32=νf_31)
    double precision :: poissonf2_12      !炭素繊維ポアソン比
    double precision :: Gf2_31            !炭素繊維剛性率(Gf_23=Gf_31)
    double precision :: nonl2
    double precision :: hf2_p1
    double precision :: hf2_p2
    double precision :: hf2_p3
    
    !時間依存ステップ
    integer swich_end
    integer step
    double precision ti
    double precision strain_giv
    
    !Von Mises
    double precision macroS_vm                                            !巨視的なフォンミーゼス応力
    double precision macroE_vm                                            !巨視的なフォンミーゼスひずみ
    double precision, allocatable, dimension(:) :: microS_vm_e       !要素ごとの微視的なフォンミーゼス応力
    double precision, allocatable, dimension(:,:) :: microS_vm_g     !積分点ごとの微視的なフォンミーゼス応力
    double precision, allocatable, dimension(:) :: microS_vm_n       !節点ごとの微視的なフォンミーゼス応力
    
    !各種マトリックス
    double precision, allocatable, dimension(:,:,:) :: Dmatrix             !各要素のDmatrix
    double precision, allocatable, dimension(:,:,:) :: Smatrix             !各軸方向の繊維に対するSmatrix
    double precision, allocatable, dimension(:,:,:) :: invSmatrix
    double precision, allocatable, dimension(:,:,:,:) :: Bmatrix
    double precision, allocatable, dimension(:,:) :: Kmatrix
    double precision, allocatable, dimension(:,:) :: Kmat_LU
    double precision, allocatable, dimension(:,:,:) :: eKmatrix
    double precision, allocatable, dimension(:,:,:) :: eKmat_cc
    double precision, allocatable, dimension(:,:,:) :: eKmat_cn
    double precision, allocatable, dimension(:,:,:) :: eKmat_nc
    double precision, allocatable, dimension(:,:,:) :: eKmat_nn
    double precision, allocatable, dimension(:,:,:) :: inveKmat_nn
    double precision, allocatable, dimension(:,:,:) :: eKcn_inveKnn
    double precision, allocatable, dimension(:,:,:) :: eKcn_inveKnn_eKnc
    double precision, allocatable, dimension(:,:,:,:) :: Hmatrix
    double precision, allocatable, dimension(:,:,:,:) :: Jmatrix
    double precision, allocatable, dimension(:,:,:,:) :: invJmatrix
    double precision, allocatable, dimension(:,:,:,:) :: cartesian
    double precision, allocatable, dimension(:,:,:,:) :: tBmatDmat
    double precision, allocatable, dimension(:,:,:,:) :: tBmatDmatBmat
    double precision, allocatable, dimension(:,:,:) :: eFmatrix
    double precision, allocatable, dimension(:,:,:) :: eFmat_c
    double precision, allocatable, dimension(:,:,:) :: eFmat_n
    double precision, allocatable, dimension(:,:,:) :: eKcn_inveKnn_eFn
    double precision, allocatable, dimension(:,:) :: Fmatrix            
    double precision, allocatable, dimension(:,:) :: Xmatrix
    double precision, allocatable, dimension(:,:,:) :: eXmat_c   
    double precision, allocatable, dimension(:,:,:) :: eKnc_eXc
    double precision, allocatable, dimension(:,:,:) :: eXmat_n
    double precision, allocatable, dimension(:,:,:) :: eXmatrix
    double precision, allocatable, dimension(:,:,:,:) :: BmatXmat
    double precision, allocatable, dimension(:,:,:,:) :: DmatBmatXmat
    double precision, allocatable, dimension(:,:,:,:) :: eAmatrix
    double precision, allocatable, dimension(:,:) :: Amatrix
    double precision, allocatable, dimension(:,:) :: invAmatrix
    
    double precision, allocatable, dimension(:,:) :: Imatrix   !単位行列
    double precision, allocatable, dimension(:,:) :: temp      !行列を一時的に保存するために使う
    double precision, allocatable, dimension(:,:) :: invtemp   !逆行列を一時的に保存するために使う
    

    double precision, allocatable, dimension(:) :: macroS              !巨視的応力
    double precision, allocatable, dimension(:) :: macroE              !巨視的ひずみ
    double precision, allocatable, dimension(:,:) :: macroE_tensor     !巨視的ひずみテンソル
    double precision, allocatable, dimension(:,:) :: microS_element    !微視的応力(要素ごと)
    double precision, allocatable, dimension(:,:,:) :: microS_gauss    !微視的応力(積分点ごと)
    double precision, allocatable, dimension(:,:) :: microS_node       !微視的応力(節点ごと)
    double precision, allocatable, dimension(:,:) :: u_sharp           !擾乱変位
    double precision, allocatable, dimension(:,:) :: microU            !微視的変位
    integer, allocatable, dimension(:) :: n_counter                      !重なっている節点の数を数える
    
    double precision, allocatable, dimension(:,:) :: hf
    double precision, allocatable, dimension(:,:,:) :: S_dev
    double precision, allocatable, dimension(:,:) :: S_equ
    double precision, allocatable, dimension(:,:,:) :: Beta
    double precision, allocatable, dimension(:,:) :: Ep_acc
    double precision, allocatable, dimension(:,:) :: eGvector
    double precision, allocatable, dimension(:,:) :: eGvec_c
    double precision, allocatable, dimension(:,:) :: eGvec_n
    double precision, allocatable, dimension(:,:) :: eKcn_inveKnn_eGn
    double precision, allocatable, dimension(:) :: Gvector
    double precision, allocatable, dimension(:) :: Pvector
    double precision, allocatable, dimension(:,:) :: ePvec_c
    double precision, allocatable, dimension(:,:) :: eKnc_ePc
    double precision, allocatable, dimension(:,:) :: ePvec_n
    double precision, allocatable, dimension(:,:) :: ePvector
    double precision, allocatable, dimension(:,:,:) :: DmatBeta
    double precision, allocatable, dimension(:,:,:,:) ::  DmatBmat
    double precision, allocatable, dimension(:,:,:) :: eRvector
    double precision, allocatable, dimension(:) :: Rvector
    double precision, allocatable, dimension(:) :: macroS_velocity
    double precision, allocatable, dimension(:) :: macroE_velocity
    double precision, allocatable, dimension(:,:,:) :: microS_velocity
    double precision, allocatable, dimension(:,:) :: u_sharp_velocity
    double precision, allocatable, dimension(:,:) :: node_after
     double precision, allocatable, dimension(:) :: Ep_el
    !LAPACK用
    character :: trans = "N"
    integer info
    integer, allocatable, dimension(:) :: ipiv
    integer, allocatable, dimension(:) :: ipiv_k !Pvectorの計算専用
    
end module