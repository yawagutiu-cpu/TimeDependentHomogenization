subroutine output_paraview

use set_parameter
implicit none

integer i,j
integer, allocatable :: offsets_arr(:)
integer, allocatable :: types_arr(:)

! offsets_arr��types_arr��ecount�̃T�C�Y�Ŋm��
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

  ! --- �C���ӏ�: offsets�̏o�� ---
  write(1,"(a)")'<DataArray type="Int32" Name="offsets" format="ascii">'
  do i=1,ecount
    offsets_arr(i) = 8 * i ! �e�I�t�Z�b�g�l��z��Ɋi�[
  enddo
  write(1,"(i0,:,1x)") offsets_arr ! �z��S�̂�1�s�ŏo��
  write(1,"(a)")'</DataArray>'

  ! --- �C���ӏ�: types�̏o�� ---
  write(1,"(a)")'<DataArray type="UInt8" Name="types" format="ascii">'
  do i=1,ecount
    types_arr(i) = 12 ! �e�^�C�v�l��z��Ɋi�[
  enddo
  write(1,"(i0,:,1x)") types_arr ! �z��S�̂�1�s�ŏo��
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

! �m�ۂ��������������
deallocate(offsets_arr)
deallocate(types_arr)

end subroutine