module mod_prkn
        implicit none
        ! Stomer-Verlet methods for testing purpose
        integer, parameter :: rkn_const_2_1a_s = 1
        real(8), parameter :: rkn_const_2_1a_aba(3) = [0.5d0, 1d0, 0.5d0]
        integer, parameter :: rkn_const_2_1_s = 1
        real(8), parameter :: rkn_const_2_1_bab(3) = [0.5d0, 1d0, 0.5d0]
        
        ! S. Blanes, P.C. Moan, ``Practical symplectic partitioned Runge-Kutta and 
        ! Runge-Kutta-Nystr\"{o}m methods,'' JCAM, 142(2002) 313-330.
        integer, parameter :: rkn_const_4_6_s = 6
        real(8), parameter :: rkn_const_4_6_bab(13) = [&
        0.0829844064174052d0,& ! cb(1)
        0.245298957184271d0,&  ! ca(1)
        0.396309801498368d0,&  ! cb(2)
        0.604872665711080d0,&  ! ca(2)
        -0.0390563049223486d0,&! cb(3)
        -0.350171622895351d0,& ! ca(3) = 0.5-sum(ca(1:2))
        0.1195241940131508d0,& ! = 1-2*sum(cb(1:3))
        -0.350171622895351d0,& ! ca(3)
        -0.0390563049223486d0,&! cb(3)
        0.604872665711080d0,&  ! ca(2)
        0.396309801498368d0,&  ! cb(2)
        0.245298957184271d0,&  ! ca(1)
        0.0829844064174052d0]  ! cb(1)
        
        integer, parameter :: rkn_const_6_11_s = 11
        real(8), parameter :: rkn_const_6_11_bab(23) = [&
        0.0414649985182624d0,&
        0.123229775946271d0,&
        0.198128671918067d0,&
        0.290553797799558d0,&
        -0.0400061921041533d0,&
        -0.127049212625417d0,&
        0.0752539843015807d0,&
        -0.246331761062075d0,&
        -0.0115113874206879d0,&
        0.357208872795928d0,&
        0.2366699247869311d0,&
        0.20477705429147d0,&
        0.2366699247869311d0,&
        0.357208872795928d0,&
        -0.0115113874206879d0,&
        -0.246331761062075d0,&
        0.0752539843015807d0,&
        -0.127049212625417d0,&
        -0.0400061921041533d0,&
        0.290553797799558d0,&
        0.198128671918067d0,&
        0.123229775946271d0,&
        0.0414649985182624d0&
        ]
        
        integer, parameter :: prk_const_6_11_s = 10
        real(8), parameter :: prk_const_6_11_bab(21) = [&
        0.0502627644003922d0,&
        0.148816447901042d0,&
        0.413514300428344d0,&
        -0.132385865767784d0,&
        0.0450798897943977d0,&
        0.067307604692185d0,&
        -0.188054853819569d0,&
        0.432666402578175d0,&
        0.541960678450780d0,&
        -0.016404589403618d0,&
        -0.7255255585086898d0,&
        -0.016404589403618d0,&
        0.541960678450780d0,&
        0.432666402578175d0,&
        -0.188054853819569d0,&
        0.067307604692185d0,&
        0.0450798897943977d0,&
        -0.132385865767784d0,&
        0.413514300428344d0,&
        0.148816447901042d0,&
        0.0502627644003922d0&
        ]
        contains
        
        
        subroutine hp_int_init2(dz, ns, nv, v, vi)
                implicit none
                integer, intent(in) :: ns, nv
                real(8), intent(in) :: dz
                complex(8), intent(inout), allocatable, dimension(:,:) :: v, vi
                call hp_int_init(dz, ns, nv, v, vi, rkn_const_2_1_s, rkn_const_2_1_bab)
        end subroutine
        
        subroutine hp_int_init4(dz, ns, nv, v, vi)
                implicit none
                integer, intent(in) :: ns, nv
                real(8), intent(in) :: dz
                complex(8), intent(inout), allocatable, dimension(:,:) :: v, vi
                call hp_int_init(dz, ns, nv, v, vi, rkn_const_4_6_s, rkn_const_4_6_bab)
        end subroutine
        
        subroutine hp_int_init6(dz, ns, nv, v, vi)
                implicit none
                integer, intent(in) :: ns, nv
                real(8), intent(in) :: dz
                complex(8), intent(inout), allocatable, dimension(:,:) :: v, vi
                call hp_int_init(dz, ns, nv, v, vi, rkn_const_6_11_s, rkn_const_6_11_bab)
        end subroutine
        
        subroutine hp_int_initd6(dz, ns, nv, v, vi)
                implicit none
                integer, intent(in) :: ns, nv
                real(8), intent(in) :: dz
                complex(8), intent(inout), allocatable, dimension(:,:) :: v, vi
                call hp_int_init(dz, ns, nv, v, vi, prk_const_6_11_s, prk_const_6_11_bab)
        end subroutine
        
        subroutine hp_int2(dz, ns, nb, nv, f, g, iv, v, vi)
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns*rkn_const_2_1_s+1), vi(nv, ns*rkn_const_2_1_s+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                call hp_int_spa(dz, ns, nb, nv, f, g, iv, v, vi, rkn_const_2_1_s, rkn_const_2_1_bab)
        end subroutine
        
        subroutine hp_int4(dz, ns, nb, nv, f, g, iv, v, vi)
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns*rkn_const_4_6_s+1), vi(nv, ns*rkn_const_4_6_s+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                call hp_int_spa(dz, ns, nb, nv, f, g, iv, v, vi, rkn_const_4_6_s, rkn_const_4_6_bab)
        end subroutine
        
        subroutine hp_int6(dz, ns, nb, nv, f, g, iv, v, vi)
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns*rkn_const_6_11_s+1), vi(nv, ns*rkn_const_6_11_s+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                call hp_int_spa(dz, ns, nb, nv, f, g, iv, v, vi, rkn_const_6_11_s, rkn_const_6_11_bab)
        end subroutine
        
        subroutine hp_int6_1d(dz, ns, nb, nv, f, g, iv, v, vi)
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns*rkn_const_6_11_s+1), vi(nv, ns*rkn_const_6_11_s+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                call hp_int_spa1d(dz, ns, f(1,1), g(1), v, vi, rkn_const_6_11_s, rkn_const_6_11_bab)
        end subroutine
        
        subroutine hp_intd2(dz, ns, nb, nv, f, g, iv, v, vi)
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns*rkn_const_2_1_s+1), vi(nv, ns*rkn_const_2_1_s+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                call hp_intd(dz, ns, nb, nv, f, g, iv, v, vi, rkn_const_2_1_s, rkn_const_2_1_bab)
        end subroutine
        
        subroutine hp_intd6(dz, ns, nb, nv, f, g, iv, v, vi)
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns*prk_const_6_11_s+1), vi(nv, ns*prk_const_6_11_s+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                call hp_intd(dz, ns, nb, nv, f, g, iv, v, vi, prk_const_6_11_s, prk_const_6_11_bab)
        end subroutine
        
        
        subroutine hp_int_init(dz, ns, nv, v, vi, s, bab)
                use mod_scpot3, only : comp_scpot
                implicit none
                integer, intent(in) :: ns, nv, s
                real(8), intent(in) :: dz, bab(2*s+1)
                complex(8), intent(inout), allocatable, dimension(:,:) :: v, vi
                
                integer :: l, k
                real(8) :: iz
                
                allocate(v(nv,s*ns+1))
                allocate(vi(nv,s*ns+1))
                
                !$omp parallel do private(iz,l,k)
                do l=1,ns
                        iz = dz*(l-1)
                        do k=1,s
                                call comp_scpot(iz,v(1,s*(l-1)+k),vi(1,s*(l-1)+k))
                                iz = iz + dz * bab(2*k)
                        end do
                end do
                iz = dz*ns
                call comp_scpot(iz,v(1,s*ns+1),vi(1,s*ns+1))
        end subroutine
        
        subroutine hp_int_init_aba(dz, ns, nv, v, vi, s, aba)
                use mod_scpot3, only : comp_scpot
                implicit none
                integer, intent(in) :: ns, nv, s
                real(8), intent(in) :: dz, aba(2*s+1)
                complex(8), intent(inout), allocatable, dimension(:,:) :: v, vi
                
                integer :: l, k
                real(8) :: iz
                
                allocate(v(nv,s*ns))
                allocate(vi(nv,s*ns))
                
                !$omp parallel do private(iz,l,k)
                do l=1,ns
                        iz = dz*(l-1)
                        do k=1,s
                                iz = iz + dz * aba(2*k-1)
                                call comp_scpot(iz,v(1,s*(l-1)+k),vi(1,s*(l-1)+k))
                        end do
                end do
        end subroutine
        
        subroutine hp_pinth(dz, nb, m, q, p, q_is_id)
                implicit none
                integer, intent(in) :: nb
                real(8), intent(in) :: dz
                complex(8), intent(in) :: m(nb, nb), q(nb, nb)
                complex(8), intent(inout) :: p(nb, nb)
                logical, intent(in) :: q_is_id

                if( q_is_id ) then
                        p = p + dz * m
                else
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0)*dz, q, nb, m, nb, (1d0,0d0), p, nb)
                end if
        end subroutine

        subroutine hp_pinth_save(dz, nb, m, q, p, psave)
                implicit none
                integer, intent(in) :: nb
                real(8), intent(in) :: dz
                complex(8), intent(in) :: m(nb, nb), q(nb, nb)
                complex(8), intent(inout) :: p(nb, nb)
                complex(8), intent(out) :: psave(nb, nb)

                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0)*dz, q, nb, m, nb, (0d0,0d0), psave, nb)
                p = p + psave
        end subroutine

        subroutine hp_pinth_usesave(nb, p, psave)
                implicit none
                integer, intent(in) :: nb
                complex(8), intent(inout) :: p(nb, nb)
                complex(8), intent(out) :: psave(nb, nb)

                p = p + psave
        end subroutine

        subroutine hp_qinth(dz, nb, q, p, q_is_id)
                implicit none
                integer, intent(in) :: nb 
                real(8), intent(in) :: dz
                complex(8), intent(in) :: p(nb, nb)
                complex(8), intent(inout) :: q(nb, nb)
                logical, intent(in) :: q_is_id

                integer :: i

                if( q_is_id ) then
                        q = dz * p
                        do i=1,nb
                                q(i,i) = q(i,i) + 1d0
                        end do
                else
                        q = q + dz * p
                end if
        end subroutine

        subroutine hp_int_spa(dz, ns, nb, nv, f, g, iv, v, vi, s, bab)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb), s
                real(8), intent(in) :: dz, bab(2*s+1)
                complex(8), intent(in) :: g(nb), v(nv, ns*s+1), vi(nv, ns*s+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8) :: t1(nb, nb+nb), psave(nb, nb)
                complex(8) :: m(nb, nb)
                real(8) :: d, ca, cb
                integer :: k, l
                logical :: x_is_id
                
                ! transform variables to reconstruct the simplest form
                call invstransh(nb, t1, f, g)
                
                ! in the rest of the process, [X^H, Y^H] is located at t1
                x_is_id = .true.
                
                do l=1,ns
                        do k=1,s
                                cb = bab(2*k-1)
                                ca = bab(2*k)
                                
                                if( k .gt. 1 .or. l .eq. 1 ) then
                                        call vvi2mate(nb, nv, v(1,s*(l-1)+k), vi(1,s*(l-1)+k), iv, g, m)
                                end if
                                if( k .eq. 1 .and. .not. x_is_id .and. l .gt. 1 ) then
                                        call hp_pinth_usesave(nb, t1(1, nb+1), psave)
                                else
                                        call hp_pinth(cb*dz, nb, m, t1, t1(1, nb+1), x_is_id)
                                end if
                                call hp_qinth(ca*dz, nb, t1, t1(1, nb+1), x_is_id)
                                x_is_id = .false.
                        end do
                        call vvi2mate(nb, nv, v(1,s*l+1), vi(1,s*l+1), iv, g, m)
                        call hp_pinth_save(bab(1)*dz, nb, m, t1, t1(1, nb+1), psave)

                        d = cond_est(nb, t1)
                        ! print *, l, d
                        if(d > 1d3 .or. l .eq. ns) then
                                call matAbdivB(nb, t1, nb)
                                x_is_id = .true.
                        else
                                x_is_id = .false.
                        end if
                end do ! l=1,ns
                
                ! back-transform of the variables
                call stransh(nb, t1, g)
                f = conjg(transpose(t1(:,(nb+1):(2*nb))))
        end subroutine
        
        subroutine hp_int_spa_aba(dz, ns, nb, nv, f, g, iv, v, vi, s, aba)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb), s
                real(8), intent(in) :: dz, aba(2*s+1)
                complex(8), intent(in) :: g(nb), v(nv, ns*s), vi(nv, ns*s)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8) :: t1(nb, nb+nb)
                complex(8) :: m(nb, nb)
                real(8) :: d, ca, cb
                integer :: k, l
                logical :: x_is_id
                
                ! transform variables to reconstruct the simplest form
                call invstransh(nb, t1, f, g)
                
                ! in the rest of the process, [X^H, Y^H] is located at t1
                
                call hp_qinth(aba(1)*dz, nb, t1, t1(1, nb+1), .true.)
                x_is_id = .false.

                do l=1,ns
                        do k=1,s
                                cb = aba(2*k)
                                if( l.ne.ns .and. k.eq.s) then
                                        ca = 2*aba(2*s+1)
                                else
                                        ca = aba(2*k+1)
                                end if

                                call vvi2mate(nb, nv, v(1,s*(l-1)+k), vi(1,s*(l-1)+k), iv, g, m)
                                call hp_pinth(cb*dz, nb, m, t1, t1(1, nb+1), x_is_id)
                                call hp_qinth(ca*dz, nb, t1, t1(1, nb+1), x_is_id)
                                x_is_id = .false.
                        end do
                        d = cond_est(nb, t1)
                        !print *, l, d
                        if(d > 1d3 .or. l .eq. ns) then
                                call matAbdivB(nb, t1, nb)
                                x_is_id = .true.
                        end if
                end do ! l=1,ns
                
                
                ! back-transform of the variables
                call stransh(nb, t1, g)
                f = conjg(transpose(t1(:,(nb+1):(2*nb))))
        end subroutine
        
        
        subroutine hp_intd(dz, ns, nb, nv, f, g, iv, v, vi, s, bab)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb), s
                real(8), intent(in) :: dz, bab(2*s+1)
                complex(8), intent(in) :: g(nb), v(nv, ns*s+1), vi(nv, ns*s+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8) :: t1(nb, nb+nb)
                complex(8) :: m(nb, nb)
                real(8) :: d, ca, cb
                integer :: k, l
                logical :: x_is_id
                
                ! transform variables to reconstruct the simplest form
                call invstransh(nb, t1, f, g)
                
                ! in the rest of the process, [X^H, Y^H] is located at t1
                x_is_id = .true.
                
                do l=1,ns
                        do k=1,s
                                if( l.ne.1 .and. k.eq.1) then
                                        cb = 2*bab(1)
                                else
                                        cb = bab(2*k-1)
                                end if
                                ca = bab(2*k)
                                
                                call vvi2matewog(nb, nv, v(1,s*(l-1)+k), vi(1,s*(l-1)+k), iv, m)
                                call hp_pinth(cb*dz, nb, m, t1, t1(1,nb+1), x_is_id)
                                call applyghi(nb, ca*dz, g, t1, t1(1,nb+1), x_is_id)
                                x_is_id = .false.
                        end do
                        d = cond_est(nb, t1)
                        !print *, l, d
                        if(d > 1d3 .or. l .eq. ns) then
                                call matAbdivB(nb, t1, nb)
                                x_is_id = .true.
                        else
                                x_is_id = .false.
                        end if
                end do ! l=1,ns
                
                call vvi2matewog(nb, nv, v(1,ns*s+1), vi(1,ns*s+1), iv, m)
                call hp_pinth(bab(1)*dz, nb, m, t1, t1(1, nb+1), .true.)
                
                ! back-transform of the variables
                call stransh(nb, t1, g)
                f = conjg(transpose(t1(:,(nb+1):(2*nb))))
        end subroutine
        
        subroutine hp_intd_aba(dz, ns, nb, nv, f, g, iv, v, vi, s, aba)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb), s
                real(8), intent(in) :: dz, aba(2*s+1)
                complex(8), intent(in) :: g(nb), v(nv, ns*s), vi(nv, ns*s)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8) :: t1(nb, nb+nb)
                complex(8) :: m(nb, nb)
                real(8) :: d, ca, cb
                integer :: k, l
                logical :: x_is_id
                
                ! transform variables to reconstruct the simplest form
                call invstransh(nb, t1, f, g)
                
                ! in the rest of the process, [X^H, Y^H] is located at t1
                
                call applyghi(nb, aba(1)*dz, g, t1, t1(1, nb+1), x_is_id)
                x_is_id = .false.
                
                do l=1,ns
                        do k=1,s
                                cb = aba(2*k)
                                if( l.ne.ns .and. k.eq.s) then
                                        ca = 2*aba(2*s+1)
                                else
                                        ca = aba(2*k+1)
                                end if
                                
                                call vvi2matewog(nb, nv, v(1,s*(l-1)+k), vi(1,s*(l-1)+k), iv, m)
                                call hp_pinth(cb*dz, nb, m, t1, t1(1, nb+1), x_is_id)
                                call applyghi(nb, ca*dz, g, t1, t1(1, nb+1), x_is_id)
                                x_is_id = .false.
                        end do
                        d = cond_est(nb, t1)
                        !print *, l, d
                        if(d > 1d3 .or. l .eq. ns) then
                                call matAbdivB(nb, t1, nb)
                                x_is_id = .true.
                        end if
                end do ! l=1,ns
                
                
                ! back-transform of the variables
                call stransh(nb, t1, g)
                f = conjg(transpose(t1(:,(nb+1):(2*nb))))
        end subroutine

        subroutine hp_int_spa1d(dz, ns, f, g, v, vi, s, bab)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, s
                real(8), intent(in) :: dz, bab(2*s+1)
                complex(8), intent(in) :: g, v(1, ns*s+1), vi(1, ns*s+1)
                complex(8), intent(inout) :: f
                
                complex(8) :: q, p, m
                real(8) :: ca, cb
                integer :: k, l
                
                ! transform variables to reconstruct the simplest form
                q = 1d0
                p = (0d0,1d0) * (1d0-conjg(f))/(1d0+conjg(f)) * conjg(g)
                
                do l=1,ns
                        do k=1,s
                                if( l.ne.1 .and. k.eq.1) then
                                        cb = 2*bab(1)
                                else
                                        cb = bab(2*k-1)
                                end if
                                ca = bab(2*k)
                                
                                m = -conjg(v(1, s*(l-1)+k) + vi(1, s*(l-1)+k)) - conjg(g*g)

                                p = p + cb * dz * m * q
                                q = q + ca * dz * p
                        end do
                        p = p / q
                        q = 1d0
                end do ! l=1,ns
                
                m = -conjg(v(1, ns*s+1) + vi(1, ns*s+1)) - conjg(g*g)
                p = p + bab(1) * dz * m * q
                
                ! back-transform of the variables
                f = (conjg(g) + (0d0,1d0)*p)/(conjg(g) + (0d0,-1d0)*p)
        end subroutine
end module
