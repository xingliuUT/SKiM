program esSkim

        use file_utils

        implicit none

        real (kind = 8), parameter :: pi = 3.141592653589793
        complex (kind = 8), parameter :: zi = (0., 1.)
        integer, parameter :: ION = 1, ELECTRON = 2, IMPURITY = 3, TRAPPEDE = 4

        ! prefix n: total number of grid points
        integer (kind = 8) :: nkz, nky, nkperp
        ! prefix i: iterator
        integer (kind = 8) :: ikz, iky, ikperp
        ! prefix d: grid spacing
        real (kind = 8) :: dkz, dky, dkperp
        ! suffix 0: offset value of a grid
        real (kind = 8) :: kz0, ky0, kperp0
        real (kind = 8), dimension(:), allocatable :: kz_grid, ky_grid, kperp_grid
        ! suffix 1: starting point of a scan
        real (kind = 8) :: kz1, ky1, kperp1
        ! index of starting point of a scan within a grid
        integer(kind = 8), dimension(:), allocatable :: ikperp_ref, ikz_ref, iky_ref

        ! TODO: add knobs for t'i = t'e, or f'z != f'i
        integer (kind = 8) :: nfprime, ntprime_i, ntprime_e, nomd
        integer (kind = 8) :: ifprime, itprime_i, itprime_e, iomd
        real (kind = 8) :: dfprime, dtprime_i, dtprime_e, domd
        real (kind = 8) :: tprime0_i, fprime0, tprime0_e, omd0
        real (kind = 8), dimension(:), allocatable :: fprime_grid, tprime_i_grid, &
                                                    & tprime_e_grid, omd_grid

        real (kind = 8) :: fprime, tprime_i, tprime_e, tprime, omd, kz, ky, kperp

        ! vi is vertical shift in complex plane of vz
        ! TODO: turn fr from an input parameter into a constant of the code
        real (kind = 8) :: vi, fr
        ! knob for each species to be adiabatic or non-adiabatic
        real (kind = 8) :: na_e, na_z, na_i
        ! lambda_debye2 non-trivial in ETG
        real (kind = 8) :: lambda_debye2
        ! Ti = Ti / Te, mu_e = me / mi (often D), mu_z = mz / mi
        real (kind = 8) :: Ti, Zeff, Z, mu_e, mu_z
        ! initial guess of mode frequency and growth rate
        real (kind = 8) :: omRe, omIm
        ! seed values for root finders
        complex (kind = 8) :: seed0, seed1, seed2
        ! spread between seed values 
        real(kind = 8) :: seedDevFrac

        ! prefix nstep: total number of steps
        ! prefix istep: iterator
        integer (kind = 8) :: nstepSec, istepSec
        integer (kind = 8) :: nstepInt, istepInt
        ! velocity integral limits
        real (kind = 8) :: vyFrom, vyTo, vzFrom, vzTo
        ! intTol is maximum error for integrals using two methods
        ! secTol is maximum error for root finders
        real (kind = 8) :: intTol, secTol

        ! omega is the variable in the GK equation
        ! root is the solution
        complex (kind = 8) :: omega, root
        ! temporary records of previously found root as seed values in scans
        !complex (kind = 8) :: root_ikperp_ref!, root_ikzx_ref!, root_ikyzx_ref
        !complex (kind = 8) :: root_ikperp_1!, root_ikzx_1!, root_ikyzx_1
        !complex (kind = 8) :: root_tmp_fl 


        ! moment for distribution function integral
        ! = 0 : density, = 1: velocity
        integer (kind = 8) :: moment = 0
        ! v_parallel
        complex (kind = 8) :: vPar
        
        ! output file units
        integer :: gammaUnit = 101, datUnit = 102

        ! trapped electrons related terms
        ! TODO: make those terms into a module
        real (kind = 8) :: omde, eps, nu0, coll_c, f_t, t_e, &
                           theta_c, theta_v, omde_avg, fr_te

        call init_file_utils
        call read_input_file
        call init_grids
        call init_output_files

        do iomd = 1, nomd
           omd = omd_grid(iomd)
           !TODO: formatted output with string/number mix
           print*, 'omd'
           print*, omd
           do ifprime = 1, nfprime
              fprime = fprime_grid(ifprime)
              print*, 'omn'
              print*, fprime
              do itprime_i = 1, ntprime_i
                 tprime_i = tprime_i_grid(itprime_i)
                 print*, 'omt_i'
                 print*, tprime_i
                    do itprime_e = 1, ntprime_e
 !                      tprime_e = tprime_e_grid(itprime_e)
                       tprime_e = tprime_i
                       print*, 'omt_e'
                       print*, tprime_e
                       call kyScan
                    end do
              end do
           end do
        end do

        call close_output_file (gammaUnit)
        call close_output_file (datUnit)
        call finish_file_utils

contains

 subroutine kyScan 

        implicit none

        complex (kind = 8), dimension(:), allocatable :: omega_ky_scan
        real (kind = 8) :: kz_gam, kz_Dmix, Dmixing_kz_max
        complex (kind = 8) :: omega_kz_max

        allocate(omega_ky_scan(nky))
        omega_ky_scan = 0. - 9999.*zi
        iky_ref = minloc(abs(ky_grid - ky1))

        do iky = iky_ref(1), nky
           ky = ky_grid(iky)

           call kzScan(kz_gam, omega_kz_max, kz_Dmix, Dmixing_kz_max)

           if (aimag(omega_kz_max) > aimag(omega_ky_scan(iky))) &
             omega_ky_scan(iky) = omega_kz_max
           !if (.not.aimag(omega_ky_scan(iky))==-9999.) &
           !  write (gammaUnit,&
           !  '(10e12.4)') tprime_i,fprime,tprime_e,omd, &
           !  kperp,ky,kz,omega_ky_scan(iky)
        end do

        do iky = iky_ref(1), 1, -1
           ky = ky_grid(iky)

           call kzScan(kz_gam, omega_kz_max, kz_Dmix, Dmixing_kz_max)

           if (aimag(omega_kz_max) > aimag(omega_ky_scan(iky))) &
             omega_ky_scan(iky) = omega_kz_max
           !if (.not.aimag(omega_ky_scan(iky))==-9999.) &
           !  write (gammaUnit,&
           !  '(10e12.4)') tprime_i,fprime,tprime_e,omd, &
           !  kperp,ky,kz,omega_ky_scan(iky)
        end do

        !root_ikyzx_ref = omega_ky_scan(iky_ref(1))
        !root_ikyzx_1 = omega_ky_scan(1)

 end subroutine

 subroutine kzScan (kz_maxgam,omega_maxgam,kz_maxDmix,Dmixing_max)

        implicit none

        real (kind = 8),intent(out) :: kz_maxgam, kz_maxDmix, Dmixing_max
        complex (kind = 8),intent(out) :: omega_maxgam
        integer (kind = 8), dimension(:), allocatable :: ikz_maxgamma
        integer (kind = 8), dimension(:), allocatable :: ikz_maxDmixing
        real (kind = 8), dimension(:), allocatable :: Dmixing_kz_scan
        complex (kind = 8), dimension(:), allocatable :: omega_kz_scan
        real (kind = 8) :: kperp_gam,kperp_Dmix,Dmixing_kperp_max
        complex (kind = 8) :: omega_kperp_max

        allocate(Dmixing_kz_scan(nkz),omega_kz_scan(nkz))
        allocate(ikz_maxgamma(1),ikz_maxDmixing(1))
        Dmixing_kz_scan = -9999.
        omega_kz_scan = 0. -9999.*zi
        ikz_ref = minloc(abs(kz_grid - kz1))

        do ikz = ikz_ref(1), nkz
           kz = kz_grid(ikz)
           print*, "ky, kz"
           print*, ky, kz

           call kperpScan(kperp_gam, omega_kperp_max, kperp_Dmix, Dmixing_kperp_max)

           if (aimag(omega_kperp_max) > aimag(omega_kz_scan(ikz))) &
              omega_kz_scan(ikz) = omega_kperp_max
           if (Dmixing_kperp_max > Dmixing_kz_scan(ikz)) &
              Dmixing_kz_scan(ikz) = Dmixing_kperp_max
           if (.not.aimag(omega_kz_scan(ikz))==-9999.) &
              write (gammaUnit,&
              '(9e12.4)') tprime_i,fprime,tprime_e,omd*ky, &
              kperp,ky,kz,omega_kz_scan(ikz)
        end do
        write (gammaUnit,*)

        do ikz = ikz_ref(1), 1, -1
           kz = kz_grid(ikz)
           print*, "ky, kz"
           print*, ky, kz

           call kperpScan(kperp_gam, omega_kperp_max, kperp_Dmix, Dmixing_kperp_max)

           if (aimag(omega_kperp_max) > aimag(omega_kz_scan(ikz))) &
              omega_kz_scan(ikz) = omega_kperp_max
           if (Dmixing_kperp_max > Dmixing_kz_scan(ikz)) &
              Dmixing_kz_scan(ikz) = Dmixing_kperp_max
           if (.not.aimag(omega_kz_scan(ikz))==-9999.) &
              write (gammaUnit,&
              '(9e12.4)') tprime_i,fprime,tprime_e,omd*ky, &
              kperp,ky,kz,omega_kz_scan(ikz)
        end do
        write (gammaUnit,*)
        
        ikz_maxgamma = maxloc(aimag(omega_kz_scan))
        ikz_maxDmixing = maxloc(Dmixing_kz_scan)
        kz_maxgam = kz_grid(ikz_maxgamma(1))
        kz_maxDmix = kz_grid(ikz_maxDmixing(1))
        omega_maxgam = omega_kz_scan(ikz_maxgamma(1))
        Dmixing_max = Dmixing_kz_scan(ikz_maxDmixing(1))

        !root_ikzx_ref = omega_kz_scan(ikz_ref(1))
        !root_ikzx_1 = omega_kz_scan(1)

 end subroutine

 subroutine kperpScan (kperp_maxgam, omega_maxgam, kperp_maxDmix, Dmixing_max)

        implicit none

        real (kind = 8), intent(out) :: kperp_maxgam, kperp_maxDmix, Dmixing_max
        complex (kind = 8), intent(out) :: omega_maxgam
        integer (kind = 8), dimension(:), allocatable :: ikperp_maxgamma
        integer (kind = 8), dimension(:), allocatable :: ikperp_maxDmixing
        real (kind = 8), dimension(:), allocatable :: Dmixing_kperp_scan
        complex (kind = 8), dimension(:), allocatable :: omega_kperp_scan
        complex (kind = 8) :: root_tmp,om_linr,om_linr_a,om_quadr,om_quadr_a
        complex (kind = 8) :: root_fl
        real (kind = 8) :: Dmixing_tmp


        allocate(Dmixing_kperp_scan(nkperp), omega_kperp_scan(nkperp))
        allocate(ikperp_maxgamma(1),ikperp_maxDmixing(1))
        Dmixing_kperp_scan = -9999.
        omega_kperp_scan = 0. -9999.*zi
        ikperp_ref = minloc(abs(kperp_grid - kperp1))
        do ikperp = 1, nkperp
            kperp = kperp_grid(ikperp)*ky
            print*, "kperp"
            print*, kperp

            call init_roots
            print*, 'init_roots'
            print*, seed0, seed1, seed2

            call rootfinder_muller(seed1,seed0,seed2,root_tmp)

            if (istepSec==nstepSec.or.isnan(real(root_tmp))) then
                print*, 'No root is found.'
                root = seed0
                omega_kperp_scan(ikperp) = 0.-9999.*zi
            else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                    isnan(aimag(root_tmp))) then
                print*, 'No root is found.'
                root = seed0
                omega_kperp_scan(ikperp) = 0.-9999.*zi
            else
                print*, "Root is found:"
                print*, root_tmp    
                omega_kperp_scan(ikperp) = root_tmp
                root = root_tmp
                call v_parallel(root_tmp,vPar)

                !if (ikperp==ikperp_ref(1).and.aimag(omega_kperp_scan(ikperp)) > aimag(root_ikperp_ref)) &
                    !root_ikperp_ref = omega_kperp_scan(ikperp)
                if (.not.aimag(omega_kperp_scan(ikz))==-9999.) &
                    write (datUnit, '(11e12.4)') tprime_i,fprime,tprime_e,&
                             ky,kz,kperp,omd,omega_kperp_scan(ikperp),&
                             vPar
            end if
        end do
        do ikperp = ikperp_ref(1), 1, -1
            kperp = kperp_grid(ikperp)*ky
            print*, "kperp"
            print*, kperp

            call init_roots
            print*, 'init_roots'
            print*, seed0, seed1, seed2

            call rootfinder_muller(seed0,seed1,seed2,root_tmp)

            if (istepSec==nstepSec.or.isnan(real(root_tmp))) then
                print*, 'No root is found.'
                root = seed0
                omega_kperp_scan(ikperp) = 0.-9999.*zi
            else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
                print*, 'No root is found.'
                root = seed0
                omega_kperp_scan(ikperp) = 0.-9999.*zi
            else
                print*, "Root is found:"
                print*, root_tmp    
                root = root_tmp
                call v_parallel(root_tmp,vPar)
                if (aimag(root_tmp)>aimag(omega_kperp_scan(ikperp))) &
                    omega_kperp_scan(ikperp) = root_tmp
                if (.not.aimag(omega_kperp_scan(ikz))==-9999.) &
                    write (datUnit, '(11e12.4)') tprime_i,fprime,tprime_e,&
                             ky,kz,kperp,omd,omega_kperp_scan(ikperp),&
                             vPar
            end if
        end do
        write (datUnit,*) 

        do ikperp = ikperp_ref(1), nkperp
            kperp = kperp_grid(ikperp)*ky
            print*, "kperp"
            print*, kperp

            call init_roots
            print*, 'init_roots'
            print*, seed0, seed1, seed2

            call rootfinder_muller(seed1,seed0,seed2,root_tmp)
            if (istepSec==nstepSec.or.isnan(real(root_tmp))) then
                print*, 'No root is found.'
                root = seed0
                omega_kperp_scan(ikperp) = 0.-9999.*zi
            else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
                print*, 'No root is found.'
                root = seed0
                omega_kperp_scan(ikperp) = 0.-9999.*zi
            else
                print*, "Root is found:"
                print*, root_tmp    
                root = root_tmp
                call v_parallel(root_tmp,vPar)
                if (aimag(root_tmp)>aimag(omega_kperp_scan(ikperp))) &
                    omega_kperp_scan(ikperp) = root_tmp
                if (.not.aimag(omega_kperp_scan(ikz))==-9999.) &
                    write (datUnit, '(11e12.4)') tprime_i,fprime,tprime_e,&
                             ky,kz,kperp,omd,omega_kperp_scan(ikperp),&
                             vPar
            end if
        end do
        write (datUnit, *)

        ikperp_maxgamma = maxloc(aimag(omega_kperp_scan))
        Dmixing_kperp_scan = aimag(omega_kperp_scan)/(kperp_grid**2+ky**2)
        ikperp_maxDmixing = maxloc(Dmixing_kperp_scan)

        kperp_maxgam = kperp_grid(ikperp_maxgamma(1))
        omega_maxgam = omega_kperp_scan(ikperp_maxgamma(1))
        kperp_maxDmix = kperp_grid(ikperp_maxDmixing(1))
        Dmixing_max = Dmixing_kperp_scan(ikperp_maxDmixing(1))

        !if (.not.aimag(omega_maxgam)==-9999.) then  
        !    root_ikperp_ref = omega_kperp_scan(ikperp_ref(1))
        !    root_ikperp_1 = omega_kperp_scan(1)
        !end if
        deallocate(Dmixing_kperp_scan,omega_kperp_scan)
        deallocate(ikperp_maxgamma,ikperp_maxDmixing)

  end subroutine

 subroutine init_roots
 
        implicit none

        complex(kind = 8) :: omega_tmp

        if (iomd==1.and.itprime_i==1.and.ifprime==1.and.iky==iky_ref(1).and.ikz==ikz_ref(1).and.ikperp==1) then
           seed0 = omRe + omIm*zi
        else 
           seed0 = root
        end if

        seed1 = seed0*(1. - seedDevFrac)
        seed2 = seed0*(1. + seedDevFrac)

 end subroutine

 subroutine v_parallel (root,vp)
        
        implicit none

        complex(kind = 8), intent(in) :: root
        complex(kind = 8), intent(out) :: vp
        integer(kind = 8) :: species
        complex(kind = 8) :: delta_ne

        print*, 'root'
        print*, root

        if (na_e==0.) then
            vp = 0.
        else
            species = ELECTRON
            moment = 1
            call vy_integral_simpson(root,vp,species,moment)
            print*, 'v_parallel'
            print*, vp, root
            moment = 0
            call vy_integral_simpson(root,delta_ne,species,moment)
            print*, 'delta ne'
            print*, 1-delta_ne
        end if

 end subroutine

 subroutine rootfinder_secant(sd1, sd2, rt)

        implicit none

        complex (kind = 8), intent(in) :: sd1, sd2
        complex (kind = 8), intent(out) :: rt
        complex (kind = 8) :: x_2, f_2, x_1, f_1
        complex (kind = 8) :: x_tmp, f_tmp

        x_1 = sd1
        call dispersion_relation(x_1, f_1)
        x_2 = sd2
        call dispersion_relation(x_2, f_2)
        print*, 'f_1, f_2'
        print*, f_1, f_2
        x_tmp = x_1
        istepSec = 0
        do while (abs(f_1) > secTol .and. abs(f_2) > secTol .and. istepSec < nstepSec) 
            x_tmp = x_1 - f_1 * ((x_1 - x_2) / (f_1 - f_2))
                call dispersion_relation(x_tmp, f_tmp)
                if (isnan(aimag(f_tmp)).or.isnan(real(f_tmp))) then
                   istepSec = nstepSec
                   exit
                end if
                print*, 'root finder'
                print*, istepSec, x_tmp, f_tmp
                f_2 = f_1
                f_1 = f_tmp
                x_2 = x_1
                x_1 = x_tmp
                istepSec = istepSec + 1
        end do
        rt = x_tmp

  end subroutine

 subroutine rootfinder_muller(sd0, sd1, sd2, rt)

        implicit none

        complex (kind = 8), intent(in) :: sd0, sd1, sd2
        complex (kind = 8), intent(out) :: rt
        complex (kind = 8) :: x_0, x_1, x_2, f_0, f_1, f_2
        complex (kind = 8) :: h_1, h_2, g_1, g_2
        complex (kind = 8) :: coeff_b, h_tmp_p, h_tmp_m
        complex (kind = 8) :: x_tmp, f_tmp, h_tmp, g_tmp

        x_0 = sd0
        call dispersion_relation(x_0, f_0)
        x_1 = sd1
        call dispersion_relation(x_1, f_1)
        x_2 = sd2
        call dispersion_relation(x_2, f_2)
        h_1 = x_1 - x_0
        h_2 = x_2 - x_1
        g_1 = (f_1 - f_0) / h_1
        !g_2 = (f_2 - f_1)/h_2
        istepSec = 0
        x_tmp = x_1
        do while (abs(f_1) > secTol.and.abs(f_2) > secTol.and.istepSec < nstepSec)
           g_2 = (f_2 - f_1) / h_2
           g_tmp = (g_2 - g_1) / (h_2 + h_1)
           coeff_b = g_2 + h_2 * g_tmp
           h_tmp_p = (-2. * f_2) / (coeff_b + sqrt(coeff_b ** 2 - 4. * g_tmp * f_2))
           h_tmp_m = (-2. * f_2) / (coeff_b - sqrt(coeff_b ** 2 - 4. * g_tmp * f_2))
           if (abs(h_tmp_p) > abs(h_tmp_m)) then
              h_tmp = h_tmp_m
           else
              h_tmp = h_tmp_p
           end if
           x_tmp = x_2 + h_tmp
           call dispersion_relation(x_tmp, f_tmp)
           if (isnan(aimag(f_tmp)) .or. isnan(real(f_tmp))) then
              istepSec = nstepSec
              exit
           end if
           x_0 = x_1
           x_1 = x_2
           x_2 = x_tmp
           f_0 = f_1
           f_1 = f_2
           f_2 = f_tmp
           h_1 = h_2
           h_2 = h_tmp
           g_1 = g_2
           istepSec = istepSec + 1
        end do
        rt = x_tmp

 end subroutine

subroutine dispersion_relation(omega, rhs)
        
        implicit none

        complex(kind = 8), intent(in) :: omega
        complex (kind = 8), intent(out) :: rhs
        complex (kind = 8) :: integral_e, integral_i, integral_z
        complex (kind = 8) :: integral_te
        integer(kind = 8) :: species

        ! 0th moment of distribution function
        moment = 0

        !trapped electron term
        if (t_e == 0.) then
            integral_te = 0.
        else
            species = TRAPPEDE
            call vy_integral_simpson(omega, integral_te, species, moment)

            print*, 'integral_te'
            print*, integral_te
        end if

        !passing electron term
        if (na_e == 0.) then
            integral_e = 0.
        else
            species = ELECTRON    !e
            call vy_integral_simpson(omega, integral_e, species, moment)
            print*, 'delta n_e'
            print*, 1 - integral_e
        end if

        !ion term
        if (na_i == 0.) then
            integral_i = 0.
        else
            species = ION    !ti
            call vy_integral_simpson(omega, integral_i, species, moment)
            print*, 'delta n_i'
            print*, 1 - integral_i
        end if

        !impurity term
        if (na_z == 0.) then
            integral_z = 0.
        else
            species = IMPURITY    !z
            call vy_integral_simpson(omega, integral_z, species, moment)
        end if

        rhs = 1. + Zeff / Ti - na_e * integral_e - &
              na_i / Ti * (Z - Zeff) / (Z - 1.) * integral_i - &
              na_z / Ti * (Zeff - 1.) * Z / (Z - 1.) * integral_z - &
              t_e * integral_te + kperp ** 2 * lambda_debye2

end subroutine

subroutine vy_integral_midpoint(omega,vy_a,vy_b,integral,species)

        implicit none

        integer(kind = 8), intent(in) :: species
        real(kind = 8),intent(in) :: vy_a, vy_b
        complex(kind = 8), intent(in) :: omega
        complex (kind = 8), intent(out) :: integral
        real(kind = 8) :: vy_m, dvy, vy_s, vy_f
        complex(kind = 8) :: f_m
        ! cut_buffer is used in dealling with damped mode
        real(kind = 8) :: cut_buffer

        cut_buffer = 0.01
        integral = 0.
        dvy = 0.01 * cut_buffer
        vy_s = vy_a
        vy_f = vy_s + dvy
        do while (vy_s < vy_b)
                vy_m = (vy_s + vy_f) / 2.
                call vz_integral(vy_m, f_m, species, omega, moment)
                integral = integral + f_m * dvy
                vy_s = vy_f
                vy_f = vy_f + dvy
                if (vy_f > vy_b) then
                    vy_f = vy_b
                    dvy = vy_f - vy_s
                end if
        end do

end subroutine

subroutine vy_integral_simpson(omega, integral, species, moment)

        implicit none

        integer(kind = 8), intent(in) :: species, moment
        complex(kind = 8), intent(in) :: omega
        complex (kind = 8), intent(out) :: integral
        real(kind = 8) :: x_left, x_center, x_right, dx, dxUpLimit
        complex (kind = 8) :: f_left, f_center, f_right
        complex (kind = 8) :: i_trapezoid, i_simpson
        integer :: vy_istep
        real(kind = 8)  :: vy_a, vy_b

        dxUpLimit = .2
        vy_a = vyFrom
        vy_b = vyTo
        vy_istep = 0.
        dx = .1
        x_left = vy_a
        x_right = x_left + dx
        call vz_integral(x_left, f_left, species, omega, moment)
        integral = 0.

        do while (x_left < vy_b .and. vy_istep < nstepInt)
                x_center = .5 * (x_left + x_right)
                call vz_integral(x_center, f_center, species, omega, moment)
                call vz_integral(x_right, f_right, species, omega, moment)
                i_trapezoid = .5 * (f_left + 2. * f_center + f_right) * .5 * dx
                i_simpson = dx * (f_left + 4. * f_center + f_right) / 6.
                do while (abs(i_trapezoid - i_simpson) > intTol)
                        dx = fr * dx * (abs(i_trapezoid - i_simpson) / intTol) ** (- 1. / 3.)
                        if (dx > dxUpLimit) dx = dxUpLimit
                        x_right = x_left + dx
                        if (x_right > vy_b) then
                                x_right = vy_b
                                dx = x_right - x_left
                        end if
                        x_center = .5 * (x_left + x_right)
                        call vz_integral(x_center, f_center, species, omega, moment)
                        call vz_integral(x_right, f_right, species, omega, moment)
                        i_trapezoid = .5 * (f_left + 2. * f_center + f_right) * .5 * dx
                        i_simpson = dx * (f_left + 4. * f_center + f_right) / 6.
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr * dx * (abs(i_trapezoid - i_simpson) /intTol) ** (- 1. / 3.)
                if (dx > dxUpLimit) dx = dxUpLimit
                x_right = x_left + dx
                if (x_right > vy_b) then
                        x_right = vy_b
                        dx = x_right - x_left
                end if
                f_left = f_right
                vy_istep = vy_istep + 1
        end do

        if (vy_istep == nstepInt) then
           istepInt = nstepInt
           integral = 0.
        end if

  end subroutine

  subroutine vz_integral(vy, integral, species, omega, moment)

        implicit none

        integer(kind = 8), intent(in) :: species, moment
        real(kind = 8), intent(in) :: vy
        complex(kind = 8), intent(in) :: omega
        real(kind = 8) :: vz_a, vz_b
        complex (kind = 8), intent(out) :: integral
        real(kind = 8) :: x_left, x_center, x_right, dx, dxUpLimit
        complex (kind = 8) :: f_left, f_center, f_right
        complex (kind = 8) :: i_trapezoid, i_simpson
        real(kind = 8) :: rho, bj, dj, fj
        integer(kind = 8) :: n
        integer(kind = 8) :: vz_istep
        real(kind = 8) :: vz_te, vy_te

        dxUpLimit = .2

        if ((.not.species == TRAPPEDE)) then
           vz_istep = 0
           vz_a = vzFrom
           vz_b = vzTo
           dx = .1
           x_left = vz_a
           x_right = x_left + dx
           f_left = integrand(x_left, vy, omega, species, moment)
           integral = 0.

           do while (x_left < vz_b .and. vz_istep < nstepInt)
                x_center = .5 * (x_left + x_right)
                f_center = integrand(x_center, vy, omega, species, moment)
                f_right = integrand(x_right, vy, omega, species, moment)
                i_trapezoid = .5 * (f_left + 2. * f_center + f_right) * .5 * dx
                i_simpson = dx * (f_left + 4. * f_center + f_right) / 6.
                do while (abs(i_trapezoid - i_simpson) > intTol)
                        dx = fr * dx * (abs(i_trapezoid - i_simpson) / intTol) ** (- 1. / 3.)
                        if (dx > dxUpLimit) dx = dxUpLimit
                        x_right = x_left + dx
                        if (x_right > vz_b) then
                                x_right = vz_b
                                dx = x_right - x_left
                        end if
                        x_center = .5 * (x_left + x_right)
                        f_center = integrand(x_center, vy, omega, species, moment)
                        f_right = integrand(x_right, vy, omega, species, moment)
                        i_trapezoid = .5 * (f_left + 2. * f_center + f_right) * .5 * dx
                        i_simpson = dx * (f_left + 4. * f_center + f_right) / 6.
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr * dx * (abs(i_trapezoid - i_simpson) / intTol) ** (- 1. / 3.)
                if (dx > dxUpLimit) dx = dxUpLimit
                x_right = x_left + dx
                if (x_right > vz_b) then 
                        x_right = vz_b
                        dx = x_right - x_left
                end if
                f_left = f_right
                vz_istep = vz_istep + 1
           end do
        else
           vy_te = vy
           vz_te = 0.
           integral = integrand(vz_te, vy_te, omega, species, moment)
        end if

        if (vz_istep == nstepInt) then
           istepInt = nstepInt
           integral = 0.
        else
           n = 0
           if (species == ION) then
              rho = kperp * vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 bj ** 2 * sqrt(1. / pi/ 2.) * vy * exp(- vy ** 2 / 2.) * integral
           else if (species == ELECTRON) then
              rho = sqrt(mu_e / Ti) * kperp * vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 bj **2 * sqrt(1. / pi / 2.) * vy * exp(- vy ** 2 / 2.) * integral
           else if (species == IMPURITY) then
              rho = sqrt(mu_z) / Z * kperp * vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 bj ** 2 *sqrt(1. / pi / 2.) * vy * exp(- vy ** 2 / 2.) * integral
           else if (species == TRAPPEDE) then
              integral = &
                 f_t * sqrt( 2. / pi) * vy ** 2 * exp(- vy ** 2 / 2.) * integral
           end if
        end if

  end subroutine

  function integrand(vz, vy, omega, species, moment)    

        real(kind = 8) :: vz, vy
        complex(kind = 8) :: omega, integrand, int_tmp
        complex(kind = 8) :: omega_star_n, omega_star_t, omega_d
        complex(kind = 8) :: omega_parallel, coll_te
        integer(kind = 8) :: species, moment

        if (species == ELECTRON .or. species == TRAPPEDE) then
           tprime = tprime_e
        else
           tprime = tprime_i
        end if

        if (species == ION) then
                omega_star_n = fprime * ky
                omega_star_t = .5 * (- 3. + vy ** 2 + &
                        vz ** 2) * tprime * ky
                omega_d = omd * ky * (.5 * vy ** 2 + vz ** 2)
                omega_parallel = kz * vz
        else if (species == ELECTRON) then
                omega_star_n = - 1. / Ti * fprime * ky
                omega_star_t = -1. / Ti * .5 * (- 3. + vy ** 2 + &
                        vz ** 2) * tprime * ky
                omega_d = - 1. / Ti * omd * ky * (.5 * vy ** 2 + vz ** 2)
                omega_parallel = sqrt(1. / Ti / mu_e) * kz * vz
        else if (species == IMPURITY) then
                omega_star_n = 1. / Z * fprime * ky
                omega_star_t = 1. / Z * .5 * (-3. + vy ** 2 + &
                        vz ** 2) * tprime * ky
                omega_d = 1. / Z * omd * ky * (.5 * vy ** 2 + vz ** 2)
                omega_parallel = sqrt(1. / mu_z) * kz * vz
        else if (species == TRAPPEDE) then
                omega_star_n = - fprime * ky / Ti
                omega_star_t = - (- 3. / 2. + vy ** 2 / 2.) &
                               * tprime * ky / Ti
                omega_d = - omde / Ti * omd * ky * vy ** 2 * omde_avg
                coll_te = coll_c * zi * 2. ** (3. /2.) * nu0 * &
                          (1. + vy / (2. + vy)) / eps / vy ** 3
        end if
        
        if (species == TRAPPEDE) then
           if (moment == 0) then
                int_tmp = (omega - omega_star_n - omega_star_t) / &
                  (omega - omega_d + coll_te)
           else if (moment == 1) then
                int_tmp = (omega - omega_star_n - omega_star_t) * vz / &
                  (omega - omega_d + coll_te)
           end if
        else
           if (moment == 0) then
                int_tmp = exp(- vz ** 2 / 2.) * &
                  (omega - omega_star_n - omega_star_t) / &
                  (omega - omega_parallel - omega_d)
           else if (moment == 1) then
                int_tmp = exp(- vz ** 2 / 2.) * vz * &
                  (omega - omega_star_n - omega_star_t) / &
                  (omega - omega_parallel - omega_d)
           end if
        end if

        integrand = int_tmp

  end function

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer(kind = 8) :: in_file
        logical :: exist

        namelist / parameters / nstepSec, nstepInt, nkz, dkz, kz0, kz1, &
                                nky, dky, ky0, ky1, dfprime, nfprime, &
                                fprime0, dtprime_i, ntprime_i, &
                                tprime0_i, dtprime_e, ntprime_e, &
                                tprime0_e, nomd, domd, omd0, &
                                vyFrom, vyTo, vzFrom, vzTo, omRe, omIm, &
                                eps, nu0, coll_c, vi, intTol, secTol, fr, &
                                Ti, Zeff, Z, mu_e, mu_z, na_e, na_z, na_i, &
                                nkperp, dkperp, kperp0, kperp1, &
                                seedDevFrac, nstepInt, t_e, omde, &
                                fr_te, lambda_debye2
        
        nstepSec = 10
        nstepInt = 1000

        fr_te = 1.
        lambda_debye2 = 0.

        nkz = 1
        dkz = .1
        kz0 = 0.
        kz1 = 0.

        nky = 1
        dky = .1
        ky0 = 0.
        ky1 = 0.
        !kperp = 1.

        dfprime = .5
        nfprime = 1.
        fprime0 = 0.

        ntprime_i = 1
        dtprime_i = 1.
        tprime0_i = 0.

        ntprime_e = 1
        dtprime_e = 1.
        tprime0_e = 0.

        nomd = 1
        domd = 1.
        omd0 = 0.

        vyFrom = 0.
        vyTo = 0.
        vzFrom = 8.0
        vzTo = 8.0

        omRe = 1.
        omIm = .1

        nu0 = 0.05
        coll_c = 1.
        eps = .2
        omde = 1.
        vi  = 1.
        secTol = 1.0E-05
        intTol = 1.0E-05
        fr = .5

        Ti = 1.
        Zeff = 1.65
        Z = 5
        mu_e = 5.45D-04
        mu_z = 10.8
        na_e = 0.
        na_z = 1.
        na_i = 1.
        t_e = 1.
        
        nkperp = 1
        dkperp = .5
        kperp0 = 0.
        kperp1 = 10.0
        
        !cut_buffer = 0..1
        seedDevFrac = 0.95
        !lower_bnd = 1E-16
        
    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit = input_unit("parameters"), nml = parameters)

    if (t_e == 0.) then
        f_t = 0.
        theta_c = pi/2.
        omde_avg = 1.
    else
        f_t = sqrt(2. * eps / (1. + eps)) * fr_te
        theta_c = acos(f_t)
        omde_avg = (3. / 4. * (pi - 2. * theta_c) + &
                   sin(2. * (pi - theta_c)) / 8. - &
                   sin(2. * theta_c) / 8.) / (pi - 2. * theta_c)
    end if

    if (t_e == 1.) then
        print*,'fraction of trapped electrons'
        print*, f_t
        print*,'critical theta'
        print*, theta_c
        print*,'average sin^2(theta)/2.+cos^2(theta)'
        print*, omde_avg
    end if

  end subroutine read_input_file

subroutine init_grids

        implicit none

        allocate(ky_grid(nky), kz_grid(nkz), kperp_grid(nkperp))

        do ikz = 1, nkz
                kz_grid(ikz) = ikz * dkz + kz0
        end do
        do iky = 1, nky
                ky_grid(iky) = iky * dky + ky0
        end do
        do ikperp = 1, nkperp
                kperp_grid(ikperp) = ikperp * dkperp + kperp0
        end do

        allocate(fprime_grid(nfprime), tprime_i_grid(ntprime_i))
        allocate(tprime_e_grid(ntprime_e))
        allocate(omd_grid(nomd))

        do ifprime = 1, nfprime
                fprime_grid(ifprime) = ifprime * dfprime + fprime0
        end do 
        do itprime_i = 1, ntprime_i
                tprime_i_grid(itprime_i) = itprime_i * dtprime_i + tprime0_i
        end do 

        do itprime_e = 1, ntprime_e
                tprime_e_grid(itprime_e) = itprime_e * dtprime_e + tprime0_e
        end do 
        do iomd = 1, nomd
                omd_grid(iomd) = iomd * domd + omd0
        end do

        allocate(ikperp_ref(1), ikz_ref(1), iky_ref(1))


end subroutine

subroutine init_output_files

        call open_output_file(gammaUnit, '.datgam')
        write (gammaUnit, '(9a12)') "tprime_i","fprime_i","tprime_e",&
                "omd","kperp","ky","kz","omega","gamma"
        call open_output_file(datUnit, '.dat')
        write (datUnit, '(11a12)') "tprime_i","fprime","tprime_e",&
                "ky","kz","kperp","omd","omega","gamma","v_parallel"
end subroutine

end program
