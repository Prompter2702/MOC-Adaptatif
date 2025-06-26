      SUBROUTINE torus_voxel(nx,ny,nz,rmaj,rmin,dx,count,voxel_grid)
         IMPLICIT NONE
         INTEGER :: nx,ny,nz
         REAL :: rmaj,rmin ! Major and minor radii
         REAL :: dx ! Voxel resolution
         INTEGER, dimension(nx,ny,nz) :: voxel_grid
         INTEGER :: i, j, k, count
         REAL :: x, y, z, dist_xy
       
         count = 0
         voxel_grid = 1
       
         do k = 1, nz
           z = (k - nz/2) * dx
           do j = 1, ny
             y = (j - ny/2) * dx
             do i = 1, nx
               x = (i - nx/2) * dx
       
               dist_xy = sqrt(x**2 + y**2)
               if (( (dist_xy - rmaj)**2 + z**2 ) <= rmin**2) then
                 voxel_grid(i,j,k) = 2
                 count = count + 1
               end if
             end do
           end do
         end do
       
         print *, 'Total filled voxels in torus:', count
         print *, 'dist_xy', count
      end SUBROUTINE torus_voxel
      
      
      
      SUBROUTINE sphere_voxel(nx,ny,nz,r,dx,voxel_grid)
         IMPLICIT NONE
         INTEGER :: nx,ny,nz
         REAL :: r! radius
         REAL :: dx ! Voxel resolution
         INTEGER, dimension(nx,ny,nz) :: voxel_grid
         INTEGER :: i, j, k, count
         REAL :: x, y, z
       
        
         count = 0
         voxel_grid = 1
        
         do k = 1, nz
           z = (k - nz/2) * dx
           do j = 1, ny
             y = (j - ny/2) * dx
             do i = 1, nx
               x = (i - nx/2) * dx
        
               if (x**2 + y**2 + z**2 <= r**2) then
                 voxel_grid(i,j,k) = 2
                 count = count + 1
               end if
             end do
           end do
         end do
        
         print *, 'Total filled voxels in sphere:', count
        
       end SUBROUTINE sphere_voxel