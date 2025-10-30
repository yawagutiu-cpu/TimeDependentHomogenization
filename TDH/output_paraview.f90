subroutine output_paraview

use set_parameter
implicit none

integer i,j
integer, allocatable :: offsets_arr(:)
integer, allocatable :: types_arr(:)

! offsets_arrとtypes_arrをecountのサイズで確保
allocate(offsets_arr(ecount))
allocate(types_arr(ecount))

open(1,file='output/result.vtu',status='replace')
  write(1,"(a)")'<?xml version="1.0"?>'
  write(1,"(a)")'<VTKFile type="UnstructuredGrid" version="1.0">'
  write(1,"(a)")'<UnstructuredGrid>'
  write(1,"(a,i0,a,i0,a)")'<Piece NumberOfPoints="', ncount, '" NumberOfCells="', ecount, '">'
  write(1,"(a)")'<Points>'
  write(1,"(a)")'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
  do i=1,ncount
    write(1,"(1p3e20.12)")(node(i,j),j=1,dim)
  enddo

  write(1,"(a)")'</DataArray>'
  write(1,"(a)")'</Points>'
  write(1,"(a)")'<Cells>'
  write(1,"(a)")'<DataArray type="Int32" Name="connectivity" format="ascii">'
  do i=1,ecount
    write(1,"(8i8)")(element(i,j)-1,j=1,ipn)
  enddo

  write(1,"(a)")'</DataArray>'

  ! --- 修正箇所: offsetsの出力 ---
  write(1,"(a)")'<DataArray type="Int32" Name="offsets" format="ascii">'
  do i=1,ecount
    offsets_arr(i) = 8 * i ! 各オフセット値を配列に格納
  enddo
  write(1,"(i0,:,1x)") offsets_arr ! 配列全体を1行で出力
  write(1,"(a)")'</DataArray>'

  ! --- 修正箇所: typesの出力 ---
  write(1,"(a)")'<DataArray type="UInt8" Name="types" format="ascii">'
  do i=1,ecount
    types_arr(i) = 12 ! 各タイプ値を配列に格納
  enddo
  write(1,"(i0,:,1x)") types_arr ! 配列全体を1行で出力
  write(1,"(a)")'</DataArray>'
  
  write(1,"(a)")'</Cells>'

  write(1,"(a)")'<PointData>'
  write(1,"(a)")'<DataArray type="Float32" Name="MicroHeni" NumberOfComponents="3" format="ascii">'
  do i=1,ncount
    write(1,"(1p3e12.4)")(microU(i,j),j=1,dim)
  enddo
  write(1,"(a)")'</DataArray>'
  write(1,"(a)")'</PointData>'
  
  write(1,"(a)")'<CellData>'
  write(1,"(a)")'<DataArray type="Float32" Name="SoutouOuryoku" NumberOfComponents="1" format="ascii">'
  do i=1,ecount
    write(1,"(1p1e12.4)") (microS_vm_e(i))
  end do
  write(1,"(a)")'</DataArray>'
  write(1,"(a)")'<DataArray type="Float32" Name="Nensose" NumberOfComponents="1" format="ascii">'
  do i=1,ecount
    write(1,"(1p1e12.4)") (real(Ep_el(i)))
  end do
  write(1,"(a)")'</DataArray>'
  write(1,"(a)")'</CellData>'

  write(1,"(a)")'</Piece>'
  write(1,"(a)")'</UnstructuredGrid>'
  write(1,"(a)")'</VTKFile>'
close(1)

! 確保したメモリを解放
deallocate(offsets_arr)
deallocate(types_arr)

end subroutine