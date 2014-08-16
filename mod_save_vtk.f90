MODULE mod_save_vtk

  CONTAINS

  subroutine save_VTK(ech,u, v,p, x, y, itsave, Nx, Ny,Psi)
  
  implicit none
  
  integer*4 Nx, Ny, i, j, itsave,ech
  real*8, dimension(0:Nx,0:Ny) :: u, v, p, Psi
  real*8, dimension(0:Nx)      :: x
  real*8, dimension(0:Ny)      :: y

  character*8 iterationname
  character*3 echantillonname
  character*40 filename
  character*256 :: foldername,full_name
  
  write(*,*)
  write(iterationname,'(I8.8)') itsave
  write(echantillonname,'(I3.3)') ech
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Création d'un nouveau dossier pour chaque échantillon
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(filename,*) echantillonname,'_champs_',iterationname,'.vtr'
  Write(foldername,*) 'Echantillon_num_',echantillonname,'.folder/'
  full_name = trim(adjustl(foldername)) // trim(adjustl(filename))
  call system('mkdir -p ' // foldername )

  open(unit=20, file=full_name, status='unknown')   
  write(20,'(A26)') '# vtk DataFile Version 2.0'
  write(20,'(A19)') 'Lid driven velocity'
  write(20,'(A5)') 'ASCII'
  write(20,'(A24)') 'DATASET RECTILINEAR_GRID'
  write(20,'((A10),3((I5)))') 'DIMENSIONS',Nx+1,Ny+1,1
  write(20,'((A13),(I5),X,(A5))') 'X_COORDINATES',Nx+1,'float'
  
  do i=0,Nx
    write(20,'(e13.6)') x(i)
  enddo
  
  write(20,'((A13),(I5),X,(A5))') 'Y_COORDINATES',Ny+1,'float'
  
  do j=0,Ny
    write(20,'(e13.6)') y(j)
  enddo
  
  write(20,'((A13),(I2),X,(A5))') 'Z_COORDINATES',1,'float'
  write(20,'(e13.6)') 0.
  write(20,'((A10),X,(I5))') 'POINT_DATA',(Nx+1)*(Ny+1)
  write(20,*) 'SCALARS Pressure float 1'
  write(20,*) 'LOOKUP_TABLE default'
  
  do j=0,Ny
    do i=0,Nx
    write(20,'(f13.6)') p(i,j)
    enddo
  enddo
  
  write(20,*) 'SCALARS Stream_function float 1'
  write(20,*) 'LOOKUP_TABLE default'
  
  do j=0,Ny
    do i=0,Nx
      write(20,'(f13.6)') Psi(i,j)
    enddo
  enddo
  
  write(20,*) 'VECTORS Velocity float'
  
  do j=0,Ny
    do i=0,Nx
      write(20,'(3(f13.6),X)') u(i,j), v(i,j), 0.
    enddo
  enddo
  
  close(20)
  
  end subroutine save_VTK

end module mod_save_vtk