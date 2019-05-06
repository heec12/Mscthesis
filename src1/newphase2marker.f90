subroutine newphase2marker (j1, j2, i1, i2, iph)
USE marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

! reset the markers within elements in the rectangular region


do kk = 1 , nmarkers
    if (mark(kk)%dead.eq.0) cycle
    n = mark(kk)%ntriag
    k = mod(n - 1, 2) + 1
    j = mod((n - k) / 2, nz-1) + 1
    i = (n - k) / 2 / (nz - 1) + 1

    if(j>=j1 .and. j<=j2 .and. i>=i1 .and. i<=i2) then
        nphase_counter(mark(kk)%phase,j,i) = nphase_counter(mark(kk)%phase,j,i) - 1
        iphase(j,i) = iph
        mark(kk)%phase = iph
        nphase_counter(iph,j,i) = nphase_counter(iph,j,i) + 1
    endif
enddo

phase_ratio(:,j,i) = 0.d0
phase_ratio(iph,j,i) = 1.d0

return
end subroutine newphase2marker


subroutine change_phase_dike
USE marker_data
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'phases.inc'
!integer ichanged(100*mnx), jchanged(100*mnx)
integer kph(1)
dimension ratio(20)

! max. depth (m) of eclogite phase transition
!change max basalt depth trying to fix the phase of partical mealting Hao Lu 5.8 2018
real*8, parameter :: max_basalt_depth = 150.e3
!real*8, parameter :: max_basalt_depth = 3.e3
! min. temperature (C) of eclogite phase transition
real*8, parameter :: min_eclogite_temp = 500.
real*8, parameter :: mantle_density = 3000.

! temperature (C) of serpentine phase transition
! change to 1550 only to see if phase fade away can be fixed (4.24.2018 Hao)
real*8, parameter :: serpentine_temp = 550.
!real*8, parameter :: serpentine_temp = 1550.
!temperature (C) and depth (m) of 10% partial melting of upper mantle.
! change to 1600 only to see if phase fade away can be fixed (4.24.2018 Hao)
real*8, parameter :: partial_melt_temp = 600.
!real*8, parameter :: partial_melt_temp = 1300.
!real*8, parameter :: partial_melt_depth = -70.e3
! thickness of new crust
real*8, parameter :: new_crust_thickness = 7.e3
!print *, new_crust_thickness
!real*8, parameter :: new_crust_thickness = 3.e3

! search the element for melting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!original code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do jj = 1, nz-1
! search for crustal depth
!dep2 = 0.25*(cord(jj,1,2)+cord(jj+1,1,2)+cord(jj,2,2)+cord(jj+1,2,2))
!if (cord(1,1,2) - dep2 >= new_crust_thickness) exit
!end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!original code ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!modified code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP parallel private(kk,i,j,k,n,tmpr,depth,iph,press,jbelow,trpres,trpres2,kinc,quad_area,yy)
!$OMP do schedule(guided)
!!!!!!!!!!!Eunseo's version!!!!!!!!!!!!
!do kk = 1 , nmarkers
 ! if (mark(kk)%dead.eq.0) cycle
  ! from ntriag, get element number
  !n = mark(kk)%ntriag
  !k = mod(n - 1, 2) + 1
  !j = mod((n - k) / 2, nz-1) + 1
  !i = (n - k) / 2 / (nz - 1) + 1

  !if( i .eq. nx/2) then
   ! if( j  .le. 6 ) then
    !  mark(kk)%phase = 3
    !else if( j .ge. 7 .and. j .lt. 15 ) then
     ! mark(kk)%phase = 5
    !else
     ! mark(kk)%phase = 4
   !end if
  !end if
!enddo
!!!!!!!!!!!!!!!Eunseo's version end!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!New version!!!!!!!!!!!!!!!!!!!!!!!!
!do ii = 1,nx-1
!do jj = 1,nz-1
!print *,'nx-1=',nx-1
!print *,'nz-1=',nz-1
!print *,'nx=',nx
!print *,'nz=',nz
!nx-1=         120
!nz-1=          40
!nx=         121 node number in x direction
!nz=          41 node number in z direction
!print *,'(nx-1)/2=',(nx)*0.5
!t = 0
do kk = 1, nmarkers !nmarkers = 43200, number of markers in the model
!print *,'(nx-1)/2=',(nx)*0.5
  !print *, 'kk =',kk
  !print *, 'nmarkers =',nmarkers
  if (mark(kk)%dead.eq.0) cycle
  !print *,'mark(kk)%dead =',mark(kk)%dead
  n = mark(kk)%ntriag
  !print *,'n = mark(kk)%ntriag =',n
  k = mod(n - 1, 2) + 1
!  print *,'mod(17,3)=', mod(17,3)
 ! print *,'mod(17,5)=', mod(17,5)
  !print *,'mod(17,3,1) +1 =', mod(17,3) + 1
  !print *,'k = mod(n - 1, 2) + 1',k
  j = mod((n - k) / 2, nz-1) + 1
  !print *,'j = mod((n - k) / 2, nz-1) + 1',j
  i = (n - k) / 2 / (nz - 1) + 1
  !print *,'time =',time
  !print *,'time_max =',time_max
  !iinj = nx/2
  !jinj = 1,nelem_inject+1 !node in z
  !print *,'i =(n - k) / 2 / (nz - 1) + 1',i
   !print *,'nx/2 =', nx/2
   !print *,'(nx-1)/2=',(nx-1)*0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!mul_phase begin!!!!!!!!!!!!!!!!!!!
temp_ave = 0.25 * (temp(j,i) + temp(j,i+1) + temp(j+1,i) + temp(j+1,i+1))
zcord_ave = 0.25 * (cord(j,i,2) + cord(j+1,i,2) + cord(j,i+1,2) + cord(j+1,i+1,2))
     if( (i .eq. (nx/2)) .and. (j .le. nelem_inject+1) ) then 
       
!       print *, ' time = ', time
!       print *, ' max_time = ', time_max

       if ((time .lt. 0.125*1e6*3.1536e7 )) then
        mark(kk)%phase = 3
       else if ((time .ge. 0.125*1e6*3.1536e7) .and. (time .lt. 0.250*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 0.250*1e6*3.1536e7) .and. (time .lt. 0.375*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 0.375*1e6*3.1536e7) .and. (time .lt. 0.500*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 0.500*1e6*3.1536e7) .and. (time .lt. 0.625*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 0.625*1e6*3.1536e7) .and. (time .lt. 0.750*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 0.750*1e6*3.1536e7) .and. (time .lt. 0.875*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 0.875*1e6*3.1536e7) .and. (time .lt. 1.0*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 1.0*1e6*3.1536e7) .and. (time .lt. 1.125*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ((time .ge. 1.125*1e6*3.1536e7) .and. (time .lt. 1.250*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 1.250*1e6*3.1536e7) .and. (time .lt. 1.375*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 1.375*1e6*3.1536e7) .and. (time .lt. 1.500*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 1.500*1e6*3.1536e7) .and. (time .lt. 1.625*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 1.625*1e6*3.1536e7) .and. (time .lt. 1.750*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 1.750*1e6*3.1536e7) .and. (time .lt. 1.875*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 1.875*1e6*3.1536e7) .and. (time .lt. 2.0*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 2.0*1e6*3.1536e7) .and. (time .lt. 2.125*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ((time .ge. 2.125*1e6*3.1536e7) .and. (time .lt. 2.250*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 2.250*1e6*3.1536e7) .and. (time .lt. 2.375*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 2.375*1e6*3.1536e7) .and. (time .lt. 2.500*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 2.500*1e6*3.1536e7) .and. (time .lt. 2.625*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 2.625*1e6*3.1536e7) .and. (time .lt. 2.750*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 2.750*1e6*3.1536e7) .and. (time .lt. 2.875*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 2.875*1e6*3.1536e7) .and. (time .lt. 3.0*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 3.0*1e6*3.1536e7) .and. (time .lt. 3.125*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ((time .ge. 3.125*1e6*3.1536e7) .and. (time .lt. 3.250*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 3.250*1e6*3.1536e7) .and. (time .lt. 3.375*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 3.375*1e6*3.1536e7) .and. (time .lt. 3.500*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 3.500*1e6*3.1536e7) .and. (time .lt. 3.625*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 3.625*1e6*3.1536e7) .and. (time .lt. 3.750*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 3.750*1e6*3.1536e7) .and. (time .lt. 3.875*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 3.875*1e6*3.1536e7) .and. (time .lt. 4.0*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 4.0*1e6*3.1536e7) .and. (time .lt. 4.125*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ((time .ge. 4.125*1e6*3.1536e7) .and. (time .lt. 4.250*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 4.250*1e6*3.1536e7) .and. (time .lt. 4.375*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 4.375*1e6*3.1536e7) .and. (time .lt. 4.500*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 4.500*1e6*3.1536e7) .and. (time .lt. 4.625*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 4.625*1e6*3.1536e7) .and. (time .lt. 4.750*1e6*3.1536e7)) then
        mark(kk)%phase = 7
       else if ((time .ge. 4.750*1e6*3.1536e7) .and. (time .lt. 4.875*1e6*3.1536e7)) then
        mark(kk)%phase = 3
       else if ( (time .ge. 4.875*1e6*3.1536e7) .and. (time .lt. 5.0*1e6*3.1536e7)) then
        mark(kk)%phase = 7

       end if
     end if

     !if( (i .eq. nx/2) .and. (j .le. nelem_inject+1) ) then
     !  if ((temp_ave.gt.600) .and. (time .lt. time_max*0.1)) then
     !   mark(kk)%phase = 5
     !  else if ((temp_ave.gt.600) .and. (time .ge. time_max*0.1) .and. (time .lt. time_max*0.2)) then
     !   mark(kk)%phase = 5
     !  else if ((temp_ave.gt.600) .and. (time .ge. time_max*0.2) .and. (time .lt. time_max*0.3)) then
     !   mark(kk)%phase = 5
     !  else if ((temp_ave.gt.600) .and. (time .ge. time_max*0.3) .and. (time .lt. time_max*0.4)) then
     !   mark(kk)%phase = 5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!mul_phase end!!!!!!!!!!!!!!!!!!!!!
end do
!$OMP end do
!$OMP end parallel


! storing plastic strain in temporary array
junk2(1:nz-1,1:nx-1) = aps(1:nz-1,1:nx-1)

! recompute phase ratio of those changed elements
do j = 1, nz-1
  i = nx/2

  kinc = sum(nphase_counter(:,j,i))
  ratio(1:nphase) = nphase_counter(1:nphase,j,i) / float(kinc)
  kph = maxloc(nphase_counter(:,j,i))

  ! the phase of this element is the most abundant marker phase
  iphase(j,i) = kph(1)
  phase_ratio(1:nphase,j,i) = ratio(1:nphase)

  ! When phase change occurs, the mineral would recrystalize and lost
  ! all plastic strain associated with this marker.
  aps(j,i) = max(aps(j,i) - junk2(j,i) / float(kinc), 0d0)
enddo

!do jj = 1, 6
!iphase(jj, nx/2) = 3
!!print *, "Basalt layer: jj=",jj," ii=", nelem_inject, " phase=", iphase(jj, nelem_inject)
!end do
!do jj = 7, 14
!iphase(jj, nx/2) = 5
!!print *, "Gabbro layer: jj=",jj," ii=", nelem_inject, " phase=", iphase(jj, nelem_inject)
!end do
!do jj = 15, nz-1
!iphase(jj, nx/2) = 4
!!print *, "Mantle layer: jj=",jj," ii=", nelem_inject, " phase=", iphase(jj, nelem_inject)
!end do

return
end subroutine change_phase_dike


