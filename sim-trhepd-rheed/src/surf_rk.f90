module mod_rk
        implicit none
        contains
        subroutine rk4_int_init(dz, ns, nv, v, vi)
                use mod_scpot3, only : comp_scpot
                implicit none
                integer, intent(in) :: ns, nv
                real(8), intent(in) :: dz
                complex(8), intent(inout), allocatable, dimension(:,:) :: v, vi
                
                integer :: l
                real(8) :: iz
                
                allocate(v(nv,2*ns+1))
                allocate(vi(nv,2*ns+1))
                
                !$omp parallel do private(iz,l)
                do l=1,ns
                        iz = dz*(l-1)
                        call comp_scpot(iz,v(1,2*(l-1)+1),vi(1,2*(l-1)+1))
                        call comp_scpot(iz+0.5d0*dz,v(1,2*(l-1)+2),vi(1,2*(l-1)+2))
                end do
                iz = dz*ns
                call comp_scpot(iz,v(1,2*ns+1),vi(1,2*ns+1))
        end subroutine
        
        subroutine rk4_int(dz, ns, nb, nv, f, g, iv, v, vi)
                use mod_matcomp
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns*2+1), vi(nv, ns*2+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8) :: t1(nb, nb+nb), t2(nb, nb+nb)
                complex(8) :: k1(nb, nb+nb), k2(nb, nb+nb), k3(nb, nb+nb), k4(nb, nb+nb)
                complex(8) :: m0(nb, nb), m1(nb, nb), m2(nb, nb)
                real(8) :: d
                integer :: l
                logical :: x_is_id
                
                ! transform variables to reconstruct the simplest form
                ! \"{y} = f(y).
                call invstransh(nb, t1, f, g)
                
                ! in the rest of the process, [X^H, Y^H] is located at t1
                call matsetid(nb, t1, nb)
                x_is_id = .true.
                
                call vvi2mate(nb, nv, v(1, 1), vi(1, 1), iv, g, m2)
                do l=1,ns
                        m0 = m2
                        call vvi2mate(nb, nv, v(1, 2*(l-1)+2), vi(1, 2*(l-1)+2), iv, g, m1)
                        call vvi2mate(nb, nv, v(1, 2*(l-1)+3), vi(1, 2*(l-1)+3), iv, g, m2)
                        k1(:,1:nb) = t1(:,(nb+1):(2*nb))
                        if( x_is_id ) then
                                k1(:,(nb+1):(2*nb)) = m0
                        else
                                call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1, nb, m0, nb, (0d0,0d0), k1(1,nb+1), nb)
                        end if
                        
                        t2 = t1 + 0.5d0 * dz * k1
                        k2(:,1:nb) = t2(:,(nb+1):(2*nb))
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t2, nb, m1, nb, (0d0,0d0), k2(1,nb+1), nb)
                        
                        t2 = t1 + 0.5d0 * dz * k2
                        k3(:,1:nb) = t2(:,(nb+1):(2*nb))
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t2, nb, m1, nb, (0d0,0d0), k3(1,nb+1), nb)
                        
                        t2 = t1 + dz * k3
                        k4(:,1:nb) = t2(:,(nb+1):(2*nb))
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t2, nb, m2, nb, (0d0,0d0), k4(1,nb+1), nb)
                        
                        t1 = t1 + (dz/6)*(k1 + 2d0*k2 + 2d0*k3 + k4)
                        x_is_id = .false.
                        
                        d = cond_est(nb, t1)
                        !print *, l, d
                        if(d > 1d3 .or. l .eq. ns) then
                                call matAbdivB(nb, t1, nb)
                                call matsetid(nb, t1, nb)
                                x_is_id = .true.
                        else
                                x_is_id = .false.
                        end if
                end do ! l=1,ns
                
                ! back-transform of the variables
                call stransh(nb, t1, g)
                f = conjg(transpose(t1(:,(nb+1):(2*nb))))
        end subroutine

        subroutine rk4_intgs(dz, ns, nb, nv, f, g, iv, v, vi)
                use mod_matcomp
                use mod_gs
                implicit none
                integer, intent(in) :: ns, nb, nv, iv(nb,nb)
                real(8), intent(in) :: dz
                complex(8), intent(in) :: g(nb), v(nv, ns*2+1), vi(nv, ns*2+1)
                complex(8), intent(inout) :: f(nb,nb)
                
                complex(8),dimension(nb, nb) :: q, p, q2, r1, s1, r2, s2, r3, s3, r4, s4 ! p2 is not used
                complex(8) :: m0(nb, nb), m1(nb, nb), m2(nb, nb), x(nb)
                real(8) :: d
                integer :: l
                
                ! transform variables to reconstruct the simplest form
                ! \"{y} = f(y).
                call invstrans(nb, q, p, f, g)
                call bcgs(nb, q, p)
                d = power_est(nb, q, p, x, -100) ! initialize x
                
                call vvi2matf(nb, nv, v(1, 1), vi(1, 1), iv, g, m2)
                do l=1,ns
                        m0 = m2
                        call vvi2matf(nb, nv, v(1, 2*(l-1)+2), vi(1, 2*(l-1)+2), iv, g, m1)
                        call vvi2matf(nb, nv, v(1, 2*(l-1)+3), vi(1, 2*(l-1)+3), iv, g, m2)

                        ! [r1; s1] = [O, I_n; m0, O] * [q; p]
                        r1 = p
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), m0, nb, q, nb, (0d0,0d0), s1, nb)
                        
                        ! [q2; p2] = [q; p] + dz/2 * [r1; s1]
                        ! [r2; s2] = [O, I_n; m1, O] * [q2; p2]
                        q2 = q + 0.5d0 * dz * r1
                        r2 = p + 0.5d0 * dz * s1
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), m1, nb, q2, nb, (0d0,0d0), s2, nb)

                        ! [q2; p2] = [q; p] + dz/2 * [r2; s2]
                        ! [r3; s3] = [O, I_n; m1, O] * [q2; p2]
                        q2 = q + 0.5d0 * dz * r2
                        r3 = p + 0.5d0 * dz * s2
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), m1, nb, q2, nb, (0d0,0d0), s3, nb)

                        ! [q2; p2] = [q; p] + dz * [r3; s3]
                        ! [r4; s4] = [O, I_n; m2, O] * [q2; p2]
                        q2 = q + dz * r3
                        r4 = p + dz * s3
                        call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), m2, nb, q2, nb, (0d0,0d0), s4, nb)
                        
                        q = q + (dz/6) * (r1 + 2d0*r2 + 2d0*r3 + r4)
                        p = p + (dz/6) * (s1 + 2d0*s2 + 2d0*s3 + s4)

                        ! using power method to estimate the spectral norm of ([q; p]^H [q; p])
                        ! which is considered to be ~ cond([q; p])
                        d = power_est(nb, q, p, x, nb)
                        if(d > 1d2) then
                                !print *, l, d
                                call bcgs(nb, q, p)
                        end if
                end do ! l=1,ns
                
                ! back-transform of the variables
                call strans(nb, q, p, g)
                f = p
        end subroutine
end module