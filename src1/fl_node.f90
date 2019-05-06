!  Calculations of forces from stresses
subroutine fl_node
use arrays

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

    !  1 - 3
    !  |   |
    !  2 - 4
    !
    !  diagonal / :
    !
    !   A:        B:
    !
    !  1---3         1
    !  | /         / |
    !  2         2---3
    !
    !  diagonal \ :
    !
    !   C:        D:
    !
    !  1          1---3
    !  | \         \  |
    !  2---3          2
    !
    !    assemblage of forces is COUNTRE CLOCK-WISE !
    !

boff = 0
! D model
Fbdyleft(:) = (/ 9.23595558e+09, -4.56360438e+09, -5.53899945e+09, -1.27145040e+10, &
    -9.01882403e+09, -2.24939998e+10, -7.29985391e+09, -1.73281850e+10, &
    -9.78016874e+09, -1.86289049e+10, -1.57499664e+10, -1.96039421e+10, &
    -3.30106903e+10, -7.82050481e+10, -1.09137224e+11, -1.07325574e+11, &
    -1.24922954e+11, -1.37466451e+11, -1.46533626e+11, -1.55405613e+11, &
    -1.63619120e+11, -1.72117329e+11, -1.80494610e+11, -1.88831597e+11, &
    -1.97133898e+11, -2.05414315e+11, -2.13680525e+11, -2.21936793e+11, &
    -2.30183708e+11, -2.38425688e+11, -2.46665916e+11, -2.54907125e+11, &
    -2.63155152e+11, -2.71415078e+11, -2.79695936e+11, -2.88000569e+11, &
    -2.96330168e+11, -3.04667706e+11, -3.12994440e+11, -3.21289757e+11, -1.62744762e+11/)
FbdyRight(:) = (/ 2.07825709e+09, -6.36432783e+09, -3.14339506e+09, -5.34381825e+09, &
    -4.01269655e+09, 5.49752421e+09, 7.58069714e+09, -1.56791921e+09, &
     6.96343206e+09, 1.13456636e+10, 2.41070892e+10, 2.96389776e+10, &
     6.01470157e+10, 8.41381458e+10, 8.47524807e+10, 1.04493067e+11, &
     1.19642993e+11, 1.39859852e+11, 1.48949368e+11, 1.56906518e+11, &
     1.64932994e+11, 1.72903639e+11, 1.81042924e+11, 1.89162347e+11, &
     1.97245555e+11, 2.05383178e+11, 2.13552232e+11, 2.21710146e+11, &
     2.29868316e+11, 2.38030646e+11, 2.46193531e+11, 2.54362321e+11, &
     2.62542466e+11, 2.70730305e+11, 2.78932170e+11, 2.87143830e+11, &
     2.95361718e+11, 3.03577195e+11, 3.11779868e+11, 3.19978203e+11, 1.62038209e+11/)

! M=0.5 model !
!FbdyLeft(:) = (/4554907782.8246479,  6608161373.5389977,  7573938504.0935974, 6842600516.6299400, &
! 4503576432.0286255, 1863402262.4682112, -981519414.23375297, -4859668044.2504082, -8536763428.7765188, &
!-12281135914.746731, -15358620876.400646, -20866223196.564064, -31853605580.192661, -31240206395.324501, &
!-31122683716.053947, -31070715154.397839, -35510684156.021690, -68348985528.507111, -113279530368.08099, &
!-141614538813.96585, -160137409886.34167, -172730502090.04483, -182834211301.44357, -191599012040.81802, &
!-199635407076.78229, -207472878611.67505, -215346891652.09903, -223281920246.69284, -231267994105.89948, &
!-239289805612.72470, -247364432998.65927, -255464502340.24866, -263595873628.99115, -271731386033.19305, &
!-279902046774.35278, -288064941304.26849, -296273034237.01666, -304411240065.67059, -312622837895.35663, &
!-320708507711.17780, -162399673526.33807/)
!FbdyRight(:) = (/1738856922.5221746, 17252201941.881931, 12243103383.099760, -4744142790.8061781, &
!-8956728547.2731533, -9056786281.7474098, -7750935970.7970991, -5452352677.7820644, -2911668732.2061510, &
!-136346910.57654953, 2783869466.9146962, 7093518645.9775677, 16824365000.435766, 19072314683.428562, &
!23228667108.316803, 26258384163.782799, 57516742289.687981,  110439612307.89200, 142568515438.55527, &
!156241775838.02243, 166820361775.56842, 175908366894.70901, 184114542351.00925, 191701884818.83551, &
!199086194235.49905, 206601370714.11414, 214340737739.47162, 222249894533.60458, 230304374646.52173, &
!238453939545.18405, 246680144442.79279, 254934554970.02621, 263188540743.69232, 271422628901.91891, &
!279698249072.17474, 287973339851.70392, 296258527342.41309, 304515637928.84808, 312791983890.00732, &
!321019075358.37268, 162332920475.08923/)

drat = dt / dt_elastic
if (drat .lt. 1.) drat = 1.

!$OMP parallel private(i, j, fx, fy, &
!$OMP                  p_est, rosubg, &
!$OMP                  press_norm_l, dlx_l, dly_l, &
!$OMP                  press_norm_r, dlx_r, dly_r, &
!$OMP                  iunknown, rho_water_g, water_depth)
!
!$OMP do
do i = 1,nx
    do j = 1,nz
        if(ynstressbc.eq.0.) then
           force(j,i,1) = 0
           force(j,i,2) = 0
           balance(j,i,1)=0
           balance(j,i,2)=0
        endif
        ! REGULAR PART - forces from stresses
        
        ! Element (j-1,i-1). Triangles B,C,D
        if ( j.ne.1 .and. i.ne.1 ) then
            ! triangle B
            ! side 2-3
            fx = stress0(j-1,i-1,1,2) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,3,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            fy = stress0(j-1,i-1,3,2) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,2,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j-1,i-1,1,2) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,3,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i-1,3,2) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,2,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle C
            ! side 2-3
            fx = stress0(j-1,i-1,1,3) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,3,3) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            fy = stress0(j-1,i-1,3,3) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,2,3) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j-1,i-1,1,3) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,3,3) * (cord(j-1,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i-1,3,3) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,2,3) * (cord(j-1,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle D
            ! side 1-2
            fx = stress0(j-1,i-1,1,4) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) - &
                 stress0(j-1,i-1,3,4) * (cord(j  ,i  ,1)-cord(j-1,i-1,1))
            fy = stress0(j-1,i-1,3,4) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) - &
                 stress0(j-1,i-1,2,4) * (cord(j  ,i  ,1)-cord(j-1,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i-1,1,4) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,3,4) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i-1,3,4) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,2,4) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif

        ! Element (j-1,i). Triangles A,B,C.
        if ( j.ne.1 .and. i.ne.nx ) then
            ! triangle A
            ! side 1-2
            fx = stress0(j-1,i  ,1,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,3,1) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            fy = stress0(j-1,i  ,3,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,2,1) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i  ,1,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,3,1) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i  ,3,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,2,1) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle B
            ! side 1-2
            fx = stress0(j-1,i  ,1,2) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 stress0(j-1,i  ,3,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            fy = stress0(j-1,i  ,3,2) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 stress0(j-1,i  ,2,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i  ,1,2) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,3,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i  ,3,2) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,2,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle C
            ! side 1-2
            fx = stress0(j-1,i  ,1,3) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,3,3) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            fy = stress0(j-1,i  ,3,3) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,2,3) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i  ,1,3) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,3,3) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i  ,3,3) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,2,3) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif
        
        ! Element (j,i-1). Triangles A,B,D
        if ( j.ne.nz .and. i.ne.1 ) then
            ! triangle A
            ! side 2-3
            fx = stress0(j  ,i-1,1,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 stress0(j  ,i-1,3,1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            fy = stress0(j  ,i-1,3,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 stress0(j  ,i-1,2,1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i-1,1,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,3,1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i-1,3,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,2,1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle B
            ! side 1-2
            fx = stress0(j  ,i-1,1,2) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,3,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i-1,3,2) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,2,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i-1,1,2) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,3,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            fy = stress0(j  ,i-1,3,2) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,2,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle D
            ! side 2-3
            fx = stress0(j  ,i-1,1,4) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,3,4) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            fy = stress0(j  ,i-1,3,4) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,2,4) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i-1,1,4) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,3,4) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i-1,3,4) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,2,4) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif

        ! Element (j,i). Triangles A,C,D
        if ( j.ne.nz .and. i.ne.nx ) then
            ! triangle A
            ! side 1-2
            fx = stress0(j  ,i  ,1,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,3,1) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i  ,3,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,2,1) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i  ,1,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,3,1) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            fy = stress0(j  ,i  ,3,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,2,1) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle C
            ! side 1-2
            fx = stress0(j  ,i  ,1,3) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,3,3) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i  ,3,3) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,2,3) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i  ,1,3) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) - &
                 stress0(j  ,i  ,3,3) * (cord(j  ,i  ,1)-cord(j+1,i+1,1))
            fy = stress0(j  ,i  ,3,3) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) - &
                 stress0(j  ,i  ,2,3) * (cord(j  ,i  ,1)-cord(j+1,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle D
            ! side 1-2
            fx = stress0(j  ,i  ,1,4) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,3,4) * (cord(j+1,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i  ,3,4) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,2,4) * (cord(j+1,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i  ,1,4) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,3,4) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            fy = stress0(j  ,i  ,3,4) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,2,4) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif

        ! GRAVITY FORCE
        force(j,i,2) = force(j,i,2) - rmass(j,i)*g
        balance(j,i,2) = balance(j,i,2) + abs(rmass(j,i)*g)

  enddo
enddo
!$OMP end do

! BOUNDARY CONDITIONS
if(nyhydro.gt.0) then
    !$OMP do
    do i=1,nx

        ! pressure from water sea on top
        rho_water_g = 1030. * g
        if(i.lt.nx) then
            water_depth = 0.5*(cord(1,i+1,2)+cord(1,i,2))
        else
            water_depth = 0.5*(cord(1,i-1,2)+cord(1,i,2))
        endif

        if (water_depth.lt.0.) then ! No water (above sea level)
            if(i.eq.1) then
                press_norm_l = 0
                dlx_l = 0
                dly_l = 0
                press_norm_r = rho_water_g*((cord(1,i+1,2)+cord(1,i,2))/2.)
                dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                dly_r = cord(1,i+1,2)-cord(1,i  ,2)
            elseif(i.eq.nx) then
                press_norm_l = rho_water_g*((cord(1,i-1,2)+cord(1,i,2))/2.)
                dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                press_norm_r = 0
                dlx_r = 0
                dly_r = 0
            else
                press_norm_l = rho_water_g*((cord(1,i-1,2)+cord(1,i,2))/2.)
                dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                press_norm_r = rho_water_g*((cord(1,i+1,2)+cord(1,i,2))/2.)
                dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                dly_r = cord(1,i+1,2)-cord(1,i  ,2)
            endif
            force(1,i,1) = force(1,i,1)-0.5*press_norm_l*dly_l-0.5*press_norm_r*dly_r
            force(1,i,2) = force(1,i,2)+0.5*press_norm_l*dlx_l+0.5*press_norm_r*dlx_r
            balance(1,i,1) = 1.0d+17
        endif
    enddo
    !$OMP end do

    !$OMP do
    do i=1,nx

        ! bottom support - Archimed force (normal to the surface, shear component = 0)
        p_est = pisos + 0.5*(den(iphsub)+drosub)*g*(cord(nz,i,2)-rzbo)
        rosubg = g * (den(iphsub)+drosub) * (1-alfa(iphsub)*temp(nz,i)+beta(iphsub)*p_est)

        if(i.eq.1) then
            press_norm_l = 0
            dlx_l = 0
            dly_l = 0

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        elseif(i.eq.nx) then
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = 0
            dlx_r = 0
            dly_r = 0
        else
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        endif
            
        force(nz,i,1) = force(nz,i,1)-0.5*press_norm_l*dly_l-0.5*press_norm_r*dly_r
        force(nz,i,2) = force(nz,i,2)+0.5*press_norm_l*dlx_l+0.5*press_norm_r*dlx_r

        balance(nz,i,1) = 1.0d+17
        !write(*,*) i,pisos,force(nz,i,1),force(nz,i,2),press_norm_l,press_norm_r,dlx_l,dlx_r,dly_l,dly_r

    enddo
    !$OMP end do
endif
!$OMP end parallel

!do j=1,nz
!    printforce(j) = force(j,1,1)
!end do
!print *, printforce

! Traction bc for side walls
!! LEFT WALL

! Fbdy calculation
!do j=1,nz
!    if ((time.lt.3.85*1.e3*3.1536e7)) then
!        SumNodeLeft(:) = 0
!        FbdyLeft(:) = 0
!    elseif ((time.gt.3.85*1.e3*3.1536e7 .and. time.lt.5*1.e3*3.1536e7)) then
!        SumNodeLeft(j) = SumNodeLeft(j) + force(j,1,1)
!    end if
!    FbdyLeft(j) = SumNodeLeft(j)/228.0
!end do
!print *, FbdyLeft

do j=1,nz
    rho_mantle_g = 1.e6 !3e8 !3300. * g

!    if ((time .lt. 40*1e3*3.1536e7)) then
!        rho_mantle_g = time*1.e6 / (1*1e6*3.1536e7)
!    if(time.lt.3.85*1e3*3.1536e7)then
!        rho_mantle_g = time*1.e8 / (1*1e6*3.1536e7)
    if ((time .lt. 11*1e3*3.1536e7)) then
        rho_mantle_g = time*1.e7 / (1*1e6*3.1536e7)

      !if () then ! For every elements in model
        if(j.eq.1) then
            press_norm_u = 0
            dlx_u = 0
            dly_u = 0
            press_norm_d = rho_mantle_g !*((cord(j+1,1,2)+cord(j,1,2))/2.)
            dlx_d = cord(j+1,1,1)-cord(j,1  ,1)
            dly_d = cord(j+1,1,2)-cord(j,1  ,2)
        elseif(j.eq.nz) then
            press_norm_u = rho_mantle_g !*((cord(j-1,1,2)+cord(j,1,2))/2.)
            dlx_u = cord(j,1  ,1)-cord(j-1,1,1)
            dly_u = cord(j,1  ,2)-cord(j-1,1,2)
            press_norm_d = 0
            dlx_d = 0
            dly_d = 0
        else
            press_norm_u = rho_mantle_g !*((cord(j-1,1,2)+cord(j,1,2))/2.)
            dlx_u = cord(j,1  ,1)-cord(j-1,1,1)
            dly_u = cord(j,1  ,2)-cord(j-1,1,2)
            press_norm_d = rho_mantle_g !*((cord(j+1,1,2)+cord(j,1,2))/2.)
            dlx_d = cord(j+1,1,1)-cord(j,1  ,1)
            dly_d = cord(j+1,1,2)-cord(j,1  ,2)
        endif
        force(j,1,1) = 0.5*press_norm_u*dly_u + 0.5*press_norm_d*dly_d
        force(j,1,2) = - 0.5*press_norm_u*dlx_u - 0.5*press_norm_d*dlx_d
        balance(j,1,1) = 1.0d+17

!    elseif ((time.gt.40*1.e3*3.1536e7 .and. temp(j,1).gt.600.0)) then
    elseif ((time.gt.3.85*1.e3*3.1536e7 .and. temp(j,1).gt.600.0)) then
        force(j,1,1) = 0.0
        force(j,1,2) = 0.0
        balance(j,1,1) = 1.0d+17

    else
        force(j,1,1) = 0.0 !force(j,1,1) - FbdyLeft(j)
        force(j,1,2) = 0.0
        balance(j,1,1) = 1.0d+17

    end if
end do

!do j=1,nz
!    printforce(j) = force(j,1,1)
!end do
!print *, printforce

!! RIGHT WALL

! Fbdy calculation
!do j=1,nz
!    if ((time.lt.3.85*1.e3*3.1536e7)) then
!        SumNodeRight(:) = 0
!        FbdyRight(:) = 0
!    elseif ((time.gt.3.85*1.e3*3.1536e7 .and. time.lt.5*1.e3*3.1536e7)) then
!        SumNodeRight(j) = SumNodeRight(j) + force(j,nx,1)
!    end if
!    FbdyRight(j) = SumNodeRight(j)/228.0
!end do

do j=1,nz
    rho_mantle_g = 1.e6 !3e8 !3300. * g

!    if ((time .lt. 40*1e3*3.1536e7)) then ! for M=0.5 f0
!        rho_mantle_g = time*1.e6 / (1*1e6*3.1536e7)
!    if ((time .lt. 3.85*1e3*3.1536e7)) then
!        rho_mantle_g = time*1.e8 / (1*1e6*3.1536e7)
    if(time.lt.11*1e3*3.1536e7)then
        rho_mantle_g = time*1.e7 / (1*1e6*3.1536e7)
        if(j.eq.1) then
            press_norm_u = 0
            dlx_u = 0
            dly_u = 0
            press_norm_d = rho_mantle_g !*((cord(j+1,nx,2)+cord(j,nx,2))/2.)
            dlx_d = cord(j+1,nx,1)-cord(j,nx  ,1)
            dly_d = cord(j+1,nx,2)-cord(j,nx  ,2)
        elseif(j.eq.nz) then
            press_norm_u = rho_mantle_g !*((cord(j-1,nx,2)+cord(j,nx,2))/2.)
            dlx_u = cord(j,nx  ,1)-cord(j-1,nx,1)
            dly_u = cord(j,nx  ,2)-cord(j-1,nx,2)
            press_norm_d = 0
            dlx_d = 0
            dly_d = 0
        else
            press_norm_u = rho_mantle_g !*((cord(j-1,nx,2)+cord(j,nx,2))/2.)
            dlx_u = cord(j,nx  ,1)-cord(j-1,nx,1)
            dly_u = cord(j,nx  ,2)-cord(j-1,nx,2)
            press_norm_d = rho_mantle_g !*((cord(j+1,nx,2)+cord(j,nx,2))/2.)
            dlx_d = cord(j+1,nx,1)-cord(j,nx  ,1)
            dly_d = cord(j+1,nx,2)-cord(j,nx  ,2)
        endif
        force(j,nx,1) = -0.5*press_norm_u*dly_u-0.5*press_norm_d*dly_d
        force(j,nx,2) = 0.5*press_norm_u*dlx_u+0.5*press_norm_d*dlx_d
        balance(j,nx,1) = 1.0d+17

!    elseif ((time.gt.40*1.e3*3.1536e7 .and. temp(j,nx).gt.600.0)) then
    elseif ((time.gt.3.85*1.e3*3.1536e7 .and. temp(j,nx).gt.600.0)) then
        force(j,nx,1) = 0.0
        force(j,nx,2) = 0.0
        balance(j,nx,1) = 1.0d+17

    else
        force(j,nx,1) = 0.0 !force(j,nx,1) - FbdyRight(j)
        force(j,nx,2) = 0.0
        balance(j,nx,1) = 1.0d+17
    end if
end do            

!do j=1,nz
!    printforce(j) = force(j,nx,1)
!end do
!print *, printforce

!$OMP do reduction(max:boff)
do i=1,nx
    do j=1,nz

        ! BALANCE-OFF
        if( iand(ncod(j,i,1),1).eq.1 .or. j.le.n_boff_cutoff ) then
            balance(j,i,1) = 0
        else
            balance(j,i,1) = abs(force(j,i,1)) / (balance(j,i,1) + 1.0d-9)
        endif

        if( iand(ncod(j,i,2),2).eq.2 .or. j.le.n_boff_cutoff ) then
            balance(j,i,2) = 0
        else
            balance(j,i,2) = abs(force(j,i,2)) / (balance(j,i,2) + 1.0d-9)
        endif

        ! DAMPING
        if( iand(ncod(j,i,1),1).ne.1 .and. abs(vel(j,i,1)).gt.1.0d-13 ) then
            force(j,i,1) = force(j,i,1) - demf*sign(force(j,i,1),vel(j,i,1))
        endif

        if( iand(ncod(j,i,2),2).ne.2 .and. abs(vel(j,i,2)).gt.1.0d-13 ) then
            force(j,i,2) = force(j,i,2) - demf*sign(force(j,i,2),vel(j,i,2))
        endif

        ! VELOCITIES FROM FORCES
        if( ncod(j,i,1) .eq. 1 ) then
            vel(j,i,1) = bc(j,i,1) 
        else
            vel(j,i,1) = vel(j,i,1) + dt*force(j,i,1)/(amass(j,i)*drat*drat)
        endif
        if( ncod(j,i,2) .eq. 1 ) then
            vel(j,i,2) = bc(j,i,2)
        else
            vel(j,i,2) = vel(j,i,2) + dt*force(j,i,2)/(amass(j,i)*drat*drat)
        endif
        ! MAX balance-off
        boff = max(boff,balance(j,i,1))
        boff = max(boff,balance(j,i,2))

    end do
end do
!$OMP end do
! Prestress to form the topo when density differences are present WITHOUT PUSHING OR PULLING!
if (i_prestress.eq.1.and.time.lt.600.e3*sec_year) then
    do k = 1,2
        do i = 1, nx
            vel(nz,i,k) = 0.
        enddo
        do j = 1, nz
            vel(j,1,k) = 0.
            vel(j,nx,k) = 0.
        enddo
    enddo
endif
return
end subroutine fl_node
