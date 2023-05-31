module mod_magnus
        implicit none
        contains
        subroutine m2_int_init(dz, ns, nv, v, vi)
                use mod_scpot3, only : comp_scpot
                implicit none
                integer, intent(in) :: ns, nv
                real(8), intent(in) :: dz
                complex(8), intent(inout), allocatable, dimension(:,:) :: v, vi
                
                integer :: l
                real(8) :: iz
                
                allocate(v(nv,ns))
                allocate(vi(nv,ns))
                
                !$omp parallel do private(iz, l)
                do l=1,ns
                        iz = dz*(l-0.5d0)
                        call comp_scpot(iz,v(1,l),vi(1,l))
                end do
        end subroutine
        
        subroutine m4_int_init(dz, ns, nv, v, vi)
                use mod_scpot3, only : comp_scpot
                implicit none
                integer, intent(in) :: ns, nv
                real(8), intent(in) :: dz
                complex(8), intent(inout), allocatable, dimension(:,:) :: v, vi
                
                integer :: l
                real(8),parameter :: c0 = (5d0-sqrt(15d0))/10, c1 = 0.5d0, c2 = (5d0+sqrt(15d0))/10
                
                allocate(v(nv,3*ns))
                allocate(vi(nv,3*ns))
                !$omp parallel do 
                do l=1,ns
                        call comp_scpot(dz*((l-1)+c0),v(1,3*l-2),vi(1,3*l-2))
                        call comp_scpot(dz*((l-1)+c1),v(1,3*l-1),vi(1,3*l-1))
                        call comp_scpot(dz*((l-1)+c2),v(1,3*l-0),vi(1,3*l-0))
                end do
        end subroutine
        
        ! Stomer-Verlet method, split XY^{-1}
        ! keep X close to identity, with threshould
        subroutine lf_int_spa(dz, ns, nb, nv, f, g, iv, v, vi)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns), vi(nv, ns)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8) :: t1(nb, nb+nb)
                complex(8) :: m(nb, nb)
                real(8) :: d, dz2
                integer :: i, l
                logical :: x_is_id
                
                ! transform variables to reconstruct the simplest form
                call invstransh(nb, t1, f, g)
                
                t1(:,1:nb) = (dz/2) * t1(:,(nb+1):(2*nb))
                do i=1,nb
                        t1(i,i) = t1(i,i) + 1d0
                end do
                x_is_id = .false.
                
                do l=1,ns
                        call vvi2mate(nb, nv, v(1,l), vi(1,l), iv, g, m)
                        
                        ! apply [I 0; xz*F I]
                        ! F <- F / (I+xz*U*F) <=>
                        ! Y^H <- Y^H + X^H*(xz*U)
                        if( x_is_id ) then
                                t1(:,(nb+1):(2*nb)) = t1(:,(nb+1):(2*nb)) + dz * m
                        else
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0)*dz, t1, nb, m, nb, (1d0,0d0), t1(1,nb+1), nb)
                        end if
                        
                        ! apply [I xz; 0 I]
                        ! F <- (xz*I+F)/I <=>
                        ! X^H <- X^H + xz*Y^H
                        if(l.eq.ns) then
                                dz2 = dz/2
                        else
                                dz2 = dz
                        endif
                        if( x_is_id ) then
                                t1(:,1:nb) = dz2 * t1(:,(nb+1):(2*nb))
                                do i=1,nb
                                        t1(i,i) = t1(i,i) + 1d0
                                end do
                        else
                                t1(:,1:nb) = t1(:, 1:nb) + dz2 * t1(:,(nb+1):(2*nb))
                        end if
                        d = cond_est(nb, t1)
                        if(d > 1d3 .or. l .eq. ns) then
                                call matAbdivB(nb, t1, nb)
                                x_is_id = .true.
                                !print *, l, d
                        else
                                x_is_id = .false.
                        end if
                end do ! l=1,ns
                
                ! back-transform of the variables
                call stransh(nb, t1, g)
                f = conjg(transpose(t1(:,(nb+1):(2*nb))))
        end subroutine
        
        ! Magnus-type 2nd-order method, split XY^{-1}
        ! keep X close to identity, with threshould
        subroutine m2_int_spa(dz, ns, nb, nv, f, g, iv, v, vi)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns), vi(nv, ns)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8) :: t1(nb, nb+nb)
                complex(8) :: m(nb, nb)
                complex(8) :: r(nb,nb), rp(nb, nb), q(nb, nb)
                real(8) :: d
                integer :: i, l
                logical :: x_is_id
                
                ! transform variables to reconstruct the simplest form
                call invstransh(nb, t1, f, g)
                
                x_is_id = .true.
                ! in the rest of the process, [X^H, Y^H] is located at t1
                
                rp = 0d0
                
                do l=1,ns
                        call vvi2mate(nb, nv, v(1,l), vi(1,l), iv, g, m)
                        
                        call hyperbolic_exp4r(dz, nb, m, nb, q, r)
                        
                        ! apply [I r+rp; 0 I]
                        ! F <- F / (I+(r+rp)*F) <=>
                        ! Y^H <- Y^H + X^H*(r^H+rp^H)
                        if( x_is_id ) then
                                t1(:,(nb+1):(2*nb)) = t1(:,(nb+1):(2*nb)) + r + rp
                        else
                                rp = rp + r
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1, nb, rp, nb, (1d0,0d0), t1(1,nb+1), nb)
                        end if
                        rp = r
                        
                        ! apply [I 0; q I]
                        ! F <- (q+F)/I <=>
                        ! X^H <- X^H + Y^H*q^H
                        if( x_is_id ) then
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1(1,nb+1), nb, q, nb, (0d0,0d0), t1, nb)
                                do i=1,nb
                                        t1(i,i) = t1(i,i) + 1d0
                                end do
                        else
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1(1,nb+1), nb, q, nb, (1d0,0d0), t1, nb)
                        end if
                        d = cond_est(nb, t1)
                        !print *, l, d
                        if(d > 1d3 .or. l .eq. ns) then
                                call matAbdivB(nb, t1, nb)
                                x_is_id = .true.
                        else
                                x_is_id = .false.
                        end if
                end do ! l=1,ns
                t1(:,(nb+1):(2*nb)) = t1(:,(nb+1):(2*nb)) + rp
                
                ! back-transform of the variables
                call stransh(nb, t1, g)
                f = conjg(transpose(t1(:,(nb+1):(2*nb))))
        end subroutine
        
        subroutine m4_int_spa(dz, ns, nb, nv, f, g, iv, v, vi)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, 3*ns), vi(nv, 3*ns)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8) :: t1(nb, nb+nb)
                complex(8) :: m1(nb, nb), m2(nb, nb), m3(nb, nb)
                complex(8) :: rp(nb, nb), r(nb,nb), q(nb, nb)
                integer :: i, l
                real(8) :: d
                logical :: x_is_id
                real(8), parameter :: c0 = sqrt(15d0)/36, c1 = 5d0/36
                
                ! transform variables to reconstruct the simplest form
                call invstransh(nb, t1, f, g)
                
                ! in the rest of the process, [X^H, Y^H] is located at t1
                
                rp = 0d0
                x_is_id = .true.
                
                do l=1,ns
                        call vvi2mate(nb, nv, v(1,3*l-2), vi(1,3*l-2), iv, g, m1)
                        call vvi2mate(nb, nv, v(1,3*l-1), vi(1,3*l-1), iv, g, m2)
                        call vvi2mate(nb, nv, v(1,3*l-0), vi(1,3*l-0), iv, g, m3)
                        r = m3 - m1
                        m3 = m1 - 2*m2 + m3
                        m1 = r
                        
                        call hyperbolic_exp4r(dz, nb, m2, nb, q, r)
                        
                        ! apply [I r+rp; 0 I]
                        ! F <- F / (I+(r+rp)*F) <=>
                        ! Y^H <- Y^H + X^H*(r^H+rp^H)
                        rp = rp + r + dz * (-c0*m1 + c1*m3)
                        if( x_is_id ) then
                                t1(:,(nb+1):(2*nb)) = t1(:,(nb+1):(2*nb)) + rp
                        else
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1, nb, rp, nb, (1d0,0d0), t1(1,nb+1), nb)
                        end if
                        
                        ! apply [I 0; q I]
                        ! F <- (q+F)/I <=>
                        ! X^H <- X^H + Y^H*q^H
                        if( x_is_id ) then
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1(1,nb+1), nb, q, nb, (0d0,0d0), t1, nb)
                                do i=1,nb
                                        t1(i,i) = t1(i,i) + 1d0
                                end do
                        else
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1(1,nb+1), nb, q, nb, (1d0,0d0), t1, nb)
                        end if
                        
                        rp = r + dz * (c0*m1 + c1*m3)
                        d = cond_est(nb, t1)
                        !print *, l, d
                        if(d > 1d3 .or. l .eq. ns) then
                                call matAbdivB(nb, t1, nb)
                                x_is_id = .true.
                        else
                                x_is_id = .false.
                        end if
                end do ! l=1,ns
                t1(:,(nb+1):(2*nb)) = t1(:,(nb+1):(2*nb)) + rp
                
                ! back-transform of the variables
                call stransh(nb, t1, g)
                f = conjg(transpose(t1(:,(nb+1):(2*nb))))
        end subroutine
        
        subroutine m6_int_spa(dz, ns, nb, nv, f, g, iv, v, vi)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, 3*ns), vi(nv, 3*ns)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8) :: t1(nb, nb+nb)
                complex(8) :: m1(nb, nb), m2(nb, nb), m3(nb, nb)
                complex(8) :: tt(nb,nb), rp(nb, nb), r(nb,nb), q(nb, nb)
                integer :: i, l
                real(8) :: d
                logical :: x_is_id
                real(8), parameter :: c0 = sqrt(15d0)/180, c1 = 1d0/18, c2 = 1d0/12960
                real(8), parameter :: d0 = 4d0/(3d0*sqrt(15d0)), d1 = 1d0/6
                
                ! transform variables to reconstruct the simplest form
                call invstransh(nb, t1, f, g)
                
                ! in the rest of the process, [X^H, Y^H] is located at t1
                
                rp = 0d0
                x_is_id = .true.
                
                do l=1,ns
                        call vvi2mate(nb, nv, v(1,3*l-2), vi(1,3*l-2), iv, g, m1)
                        call vvi2mate(nb, nv, v(1,3*l-1), vi(1,3*l-1), iv, g, m2)
                        call vvi2mate(nb, nv, v(1,3*l-0), vi(1,3*l-0), iv, g, m3)
                        r = m3 - m1
                        q = m1 - 2*m2 + m3
                        m1 = r
                        m2 = m2 + d1*q
                        
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0)*c2*(dz**3), m1, nb, m1, nb, (0d0,0d0), m3, nb)
                        m3 = m3 + c1*dz*q
                        
                        tt = m2 - d0*m1
                        call hyperbolic_exp6r(dz/2, nb, tt, q, r)
                        
                        ! apply [I r+rp; 0 I]
                        ! F <- F / (I+(r+rp)*F) <=>
                        ! Y^H <- Y^H + X^H*(r^H+rp^H)
                        rp = rp + r - (c0*dz) * m1 + m3
                        if( x_is_id ) then
                                t1(:,(nb+1):(2*nb)) = t1(:,(nb+1):(2*nb)) + rp
                        else
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1, nb, rp, nb, (1d0,0d0), t1(1,nb+1), nb)
                        end if
                        
                        ! apply [I 0; q I]
                        ! F <- (q+F)/I <=>
                        ! X^H <- X^H + Y^H*q^H
                        if( x_is_id ) then
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1(1,nb+1), nb, q, nb, (0d0,0d0), t1, nb)
                                do i=1,nb
                                        t1(i,i) = t1(i,i) + 1d0
                                end do
                        else
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1(1,nb+1), nb, q, nb, (1d0,0d0), t1, nb)
                        end if
                        
                        rp = r
                        tt = m2 + d0*m1
                        call hyperbolic_exp6r(dz/2, nb, tt, q, r)
                        
                        rp = rp + r
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1, nb, rp, nb, (1d0,0d0), t1(1,nb+1), nb)
                        
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1(1,nb+1), nb, q, nb, (1d0,0d0), t1, nb)
                        
                        rp = r + (dz*c0) * m1 + m3
                        d = cond_est(nb, t1)
                        !print *, l, d
                        if(d > 1d3 .or. l .eq. ns) then
                                call matAbdivB(nb, t1, nb)
                                x_is_id = .true.
                        else
                                x_is_id = .false.
                        end if
                end do ! l=1,ns
                t1(:,(nb+1):(2*nb)) = t1(:,(nb+1):(2*nb)) + rp
                
                ! back-transform of the variables
                call stransh(nb, t1, g)
                f = conjg(transpose(t1(:,(nb+1):(2*nb))))
        end subroutine
        
end module