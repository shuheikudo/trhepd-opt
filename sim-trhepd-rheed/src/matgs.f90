module mod_gs
        implicit none
        contains
        subroutine bcgs(n, q, p)
                implicit none
                integer, intent(in) :: n
                complex(8), intent(inout) :: q(n, n), p(n, n)
                complex(8) :: t(n+n, n)

                t(1:n, :) = q
                t((n+1):(2*n), :) = p
                call cqr(2*n, n, t)
                q = t(1:n, :)
                p = t((n+1):(2*n), :)
        end subroutine

        subroutine bcgsh(n, t)
                implicit none
                integer, intent(in) :: n
                complex(8), intent(inout) :: t(n, n+n)

                call cqrh(n, n+n, t)
        end subroutine

        subroutine cqr(m, n, a)
                implicit none
                integer, intent(in) :: m, n
                complex(8), intent(inout) :: a(m, n)

                integer :: info
                complex(8) :: c(n, n)

                call zherk('U', 'C', n, m, (1d0, 0d0), a, m, (0d0, 0d0), c, n)
                call zpotrf('U', n, c, n, info)
                ! if (info) abort
                call ztrsm('R', 'U', 'N', 'N', m, n, (1d0, 0d0), c, n, a, m)
        end subroutine

        subroutine cqrh(m, n, a)
                implicit none
                integer, intent(in) :: m, n
                complex(8), intent(inout) :: a(m, n)

                integer :: info
                complex(8) :: c(m, m)

                call zherk('U', 'N', m, n, (1d0, 0d0), a, m, (0d0, 0d0), c, m)
                call zpotrf('U', m, c, m, info)
                ! if (info) abort
                call ztrsm('L', 'U', 'C', 'N', m, n, (1d0, 0d0), c, m, a, m)
        end subroutine
end module