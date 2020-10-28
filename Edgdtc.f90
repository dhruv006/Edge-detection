program edgedetection
implicit none

integer, parameter :: dp=selected_real_kind(15,300)
integer, dimension(3,3) :: GX, GY 
integer, dimension(:,:), allocatable :: out_image, in_image
integer :: kx, ky, i, j, sx, sy, M, N, greyscale
integer :: istat
character(len=2) :: p2

open(file = "clown.pgm",unit = 11, iostat = istat)
if (istat /= 0) print *, 'error opening clown'

read(11, fmt = *)
read(11, '(i3)') M
read(11, '(i3)') N !size of image
read(11, fmt = *) greyscale 

!convolution kernels
GX = reshape((/ -1, -2, -1, 0, 0, 0, 1, 2, 1 /), (/ 3, 3 /))
GY = reshape((/ -1, 0, 1, -2, 0, 2, -1, 0, 1 /), (/ 3, 3 /))

allocate(out_image(M,N))
allocate(in_image(M,N))
!read image in
do i = 1, M
    do j=1,N
      read(11, fmt=*, iostat = istat) in_image(i,j) !pixels
      if (istat /= 0) print *, 'error reading clown.pgm'
    !print *, i, j, in_image
    end do
end do

!apply GX and GY through spatial convolution
do j = 2, N-1
    do i = 2, M-1  
        sx = 0
        sy = 0  
        do ky = 1, 3
            do kx = 1, 3
                sx = sx + in_image(i + kx - 2, j + ky - 2) * GX(kx, ky)
                sy = sy + in_image(i + kx - 2, j + ky - 2) * GY(kx, ky)
            end do
        end do
        out_image(i, j) = int(sqrt(real(sx, dp) ** 2 + real(sy, dp) ** 2)) !Find gradient
    end do
end do


!normalize
out_image = real(255*(out_image-minval(out_image)),dp)/real(maxval(out_image)-minval(out_image),dp)
close(11, iostat = istat)
if (istat /= 0) print *, 'error closing clown'
open(file = 'edgyclown.pgm', unit = 12, iostat = istat)
if (istat /= 0) print *, 'error opening edgyclown'

write(12, fmt='(a2)', iostat = istat) 'P2' 
write(12, fmt='(2(i3,1x))', iostat = istat) M, N
write(12, fmt='(i3)', iostat = istat) greyscale

!write to file
do i = 1, M
    do j = 1, N
        write (12, fmt='(i3)', iostat = istat) out_image(i,j)
        if (istat /= 0) print *, 'error writing out_image'
    end do
end do

close(12, iostat = istat)
if (istat /= 0) print *, 'error closing edgyclown.pgm'

end program edgedetection