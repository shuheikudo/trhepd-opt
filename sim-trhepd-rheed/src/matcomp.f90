module mod_matcomp
        implicit none
        contains
        
        subroutine matAdivB2(n, A, B)
                ! compute A / B
                ! where A and B are n times n matrices
                ! It will broke both A and B.
                ! A is overwriten by the result end of the subroutine.
                implicit none
                integer, intent(in) :: n
                complex(kind(0d0)), intent(inout) :: A(n, n), B(n, n)
                
                complex(kind(0d0)) :: t1(n, n+n)
                t1(:, 1:n) = transpose(B)
                t1(:, (n+1):(n+n)) = transpose(A)
                call matAbdivB(n, t1, n)
                A = transpose(t1(:, (n+1):(n+n)))
                return
        end subroutine
        
        subroutine matAdivB(n, A, lda)
                ! compute A / B
                ! where A and B are n times n matrices
                ! and B is located at A(1:n, (n+1):(2*n)).
                ! It will broke both A and B.
                ! A is overwriten by the result end of the subroutine.
                implicit none
                integer, intent(in) :: n, lda
                complex(kind(0d0)), intent(inout) :: A(lda, 2*n)
                integer :: i, j
                complex(kind(0d0)) :: t
                
                ! A*inv(B) = X  <=> B^T X^T = A^T
                ! Thus, we first transpose-swap A and B, then, call matAbdivB,
                ! and, transpose the result
                do j=1,n
                        do i=1,n
                                t = A(i, j)
                                A(i, j) = A(j, i+n)
                                A(j, i+n) = t
                        end do
                end do
                call matAbdivB(n, A, lda)
                do j=1, n
                        do i=1, n
                                A(i, j) = A(j, i+n)
                        end do
                end do
                return
        end subroutine
        
        subroutine matAbdivB(n, A, lda)
                ! compute A \ B (matlab notation, which is equivalent with
                ! inv(A)*B).
                ! where A and B are n times n matrices and B is located at 
                ! A(1:n, (n+1):(2*n)).
                ! It will broke both A and B.
                ! B is overwriten by the result end of the subroutine.
                implicit none
                integer, intent(in) :: n, lda
                complex(kind(0d0)), intent(inout) :: A(lda, 2*n)
                integer :: p(n), info
                
                !call zgesv(n, n, A, lda, p, A(1,n+1), lda, info)
                call zgetrf(n, 2*n, A, lda, p, info)
                call ztrsm('L', 'U', 'N', 'N', n, n, (1d0,0d0), A, lda, A(1,n+1), lda)
                return
        end subroutine
        
        subroutine matsetid(n, A, lda)
                implicit none
                integer, intent(in) :: n, lda
                complex(8), intent(out) :: A(n, lda)
                integer :: i, j
                !$omp parallel do private(i,j)
                do j=1,n
                        do i=1,n
                                if(i.eq.j) then
                                        A(i,j) = 1d0
                                else
                                        A(i,j) = 0d0
                                end if
                        end do
                end do
                !$omp end parallel do
        end subroutine
        
        subroutine invstrans(n, q, p, f, g)
                ! apply t1 <- S^{-1} * [I_n; f]
                ! where S = 1/2 [G, -iI_n; G, iI_n]
                ! -> S^{-1} = 2 [G^{-1}, G^{-1}; iI_n, -iI_n]
                implicit none
                integer,intent(in) :: n
                complex(8),intent(in) :: f(n, n), g(n)
                complex(8),intent(inout) :: q(n, n), p(n, n)
                
                integer :: i, j
                complex(8) :: ig
                
                !$omp parallel do private(i,j,ig)
                do j=1,n
                        do i=1,n
                                ig = 1d0 / g(i)
                                if(i .eq. j) then
                                        q(i, j) = ig + ig * f(i, j)
                                        p(i, j) = (0d0, -1d0) + (0d0, 1d0) * f(i, j)
                                else
                                        q(i, j) = ig * f(i, j)
                                        p(i, j) = (0d0, 1d0) * f(i, j)
                                end if
                        end do
                end do
                !$omp end parallel do
        end subroutine
        
        subroutine strans(n, q, p, g)
                ! apply t1 <- S * t1
                ! where S = 1/2 [\Gamma, -iI_n; \Gamma, iI_n]
                implicit none
                integer,intent(in) :: n
                complex(8),intent(in) :: g(n)
                complex(8),intent(inout) :: q(n, n), p(n, n)
                
                integer :: i, j
                complex(8) :: x, y
                !$omp parallel do private(i,j,x,y)
                do j=1,n
                        do i=1,n
                                x = q(i, j)
                                y = p(i, j)
                                q(i, j) = g(i) * x + (0d0, 1d0) * y
                                p(i, j) = g(i) * x - (0d0, 1d0) * y
                        end do
                end do
                !$omp end parallel do
                call matAdivB2(n, p, q)
        end subroutine
        
        subroutine invstransh(n, t1, f, g)
                ! Hermitian version of the above
                implicit none
                integer,intent(in) :: n
                complex(8),intent(in) :: f(n, n), g(n)
                complex(8),intent(inout) :: t1(n, n+n)
                
                integer :: i, j
                
                !$omp parallel do private(i,j)
                do j=1,n
                        do i=1,n
                                if(i .eq. j) then
                                        t1(i, j) = 1d0 + conjg(f(j, j))
                                        t1(i, n+j) = 1d0 - conjg(f(j, j))
                                else
                                        t1(i, j) = conjg(f(j, i))
                                        t1(i, n+j) = -conjg(f(j, i))
                                end if
                        end do
                end do
                !$omp end parallel do

                call matAbdivB(n, t1, n)

                !$omp parallel do private(i,j)
                do j=1,n
                        do i=1,n
                                t1(i, n+j) = (0d0, 1d0) * t1(i, n+j) * conjg(g(i))
                        end do
                end do
                !$omp end parallel do
        end subroutine
        
        subroutine stransh(n, t1, g)
                implicit none
                integer,intent(in) :: n
                complex(8),intent(in) :: g(n)
                complex(8),intent(inout) :: t1(n, n+n)
                
                integer :: i, j
                
                !$omp parallel do private(i,j)
                do j=1,n
                        do i=1,n
                                if(i .eq. j) then
                                        t1(i, j) = conjg(g(i)) + (0d0, -1d0) * t1(i, j+n)
                                        t1(i, j+n) = conjg(g(i)) + (0d0, 1d0) * t1(i, j+n)
                                else
                                        t1(i, j) = (0d0, -1d0) * t1(i, j+n)
                                        t1(i, j+n) = (0d0, 1d0) * t1(i, j+n)
                                end if
                        end do
                end do
                !$omp end parallel do
                call matAbdivB(n, t1, n)
        end subroutine
        
        
        subroutine vvi2mat(n, nv, v, vi, iv, g, t)
                implicit none
                integer,intent(in) :: n, nv, iv(n, n)
                complex(8),intent(in) :: v(nv), vi(nv), g(n)
                complex(8),intent(out) :: t(n, n)
                
                integer :: i, j
                complex(8) :: v0(nv), v1(nv)
                !$omp parallel private(i,j)
                !$omp workshare
                v0 = v + vi
                v1 = dconjg(v-vi)
                !$omp end workshare
                !$omp do
                do j=1,n
                        do i=1,j-1
                                t(i, j) = v1(iv(j,i))
                        end do
                        t(j, j) = v0(1) + g(j)*g(j)
                        do i=j+1,n
                                t(i, j) = v0(iv(i,j))
                        end do
                end do
                !$omp end do
                !$omp end parallel
        end subroutine
        
        subroutine vvi2math(n, nv, v, vi, iv, g, t)
                implicit none
                integer,intent(in) :: n, nv, iv(n, n)
                complex(8),intent(in) :: v(nv), vi(nv), g(n)
                complex(8),intent(out) :: t(n, n)
                
                integer :: i, j
                complex(8) :: v0(nv), v1(nv)
                !$omp parallel private(i,j)
                !$omp workshare
                v0 = conjg(v + vi)
                v1 = v-vi
                !$omp end workshare
                !$omp do
                do j=1,n
                        do i=1,j-1
                                t(i, j) = v0(iv(j,i))
                        end do
                        t(j, j) = v0(1) + conjg(g(j)*g(j))
                        do i=j+1,n
                                t(i, j) = v1(iv(i,j))
                        end do
                end do
                !$omp end do
                !$omp end parallel
        end subroutine
        
        subroutine vvi2mate(n, nv, v, vi, iv, g, t)
                implicit none
                integer,intent(in) :: n, nv, iv(n, n)
                complex(8),intent(in) :: v(nv), vi(nv), g(n)
                complex(8),intent(out) :: t(n, n)
                
                integer :: i, j
                complex(8) :: v0(nv), v1(nv)

                !$omp parallel private(i,j)
                !$omp workshare
                v0 = -conjg(v + vi)
                v1 = -v+vi
                !$omp end workshare
                !$omp do
                do j=1,n
                        do i=1,j-1
                                t(i, j) = v0(iv(j,i))
                        end do
                        t(j, j) = v0(1) - conjg(g(j)*g(j))
                        do i=j+1,n
                                t(i, j) = v1(iv(i,j))
                        end do
                end do
                !$omp end do
                !$omp end parallel
        end subroutine
        
        subroutine vvi2matf(n, nv, v, vi, iv, g, t)
                implicit none
                integer,intent(in) :: n, nv, iv(n, n)
                complex(8),intent(in) :: v(nv), vi(nv), g(n)
                complex(8),intent(out) :: t(n, n)
                
                integer :: i, j
                complex(8) :: v0(nv), v1(nv)
                !$omp parallel private(i,j)
                !$omp workshare
                v0 = -(v + vi)
                v1 = conjg(vi-v)
                !$omp end workshare
                !$omp do
                do j=1,n
                        do i=1,j-1
                                t(i, j) = v1(iv(j,i))
                        end do
                        t(j, j) = v0(1) - g(j)*g(j)
                        do i=j+1,n
                                t(i, j) = v0(iv(i,j))
                        end do
                end do
                !$omp end do
                !$omp end parallel
        end subroutine
        
        subroutine vvi2matewog(n, nv, v, vi, iv, t)
                ! without g
                implicit none
                integer,intent(in) :: n, nv, iv(n, n)
                complex(8),intent(in) :: v(nv), vi(nv)
                complex(8),intent(out) :: t(n, n)
                
                integer :: i, j
                complex(8) :: v0(nv), v1(nv)
                !$omp parallel private(i,j)
                !$omp workshare
                v0 = -conjg(v + vi)
                v1 = -v+vi
                !$omp end workshare
                !$omp do
                do j=1,n
                        do i=1,j-1
                                t(i, j) = v0(iv(j,i))
                        end do
                        t(j, j) = v0(1)
                        do i=j+1,n
                                t(i, j) = v1(iv(i,j))
                        end do
                end do
                !$omp end do
                !$omp end parallel
        end subroutine
        
        subroutine applyghi(n, h, g, x, y, x_is_id)
                implicit none
                integer, intent(in) :: n
                logical, intent(in) :: x_is_id
                real(8), intent(in) :: h
                complex(8), intent(in) :: g(n)
                complex(8), intent(inout) :: x(n, n), y(n, n)
                ! apply exp(h[0 G^2; I 0])^H from right
                ! exp(h[0 x^2; 1 0]) = [cosh(hx) xsinh(hx); sinh(hx)/x cosh(hx)]
                complex(8) :: d1(n), d2(n), d3(n), t0, t1
                integer :: i, j
                !$omp parallel private(i,j,t0,t1)
                !$omp do
                do i=1,n
                        t0 = h * conjg(g(i))
                        t1 = sinh(t0)
                        d1(i) = cosh(t0)
                        d2(i) = t1 * conjg(g(i))
                        d3(i) = t1 / conjg(g(i))
                end do
                !$omp end do
                if( x_is_id ) then
                        !$omp do
                        do j=1,n
                                do i=1,n
                                        if(i.eq.j) then
                                                t1 = y(i, j)
                                                y(i, j) = t1 * d1(j) + d2(j)
                                                x(i, j) = t1 * d3(j) + d1(j)
                                        else
                                                t1 = y(i, j)
                                                y(i, j) = t1 * d1(j)
                                                x(i, j) = t1 * d3(j)
                                        end if
                                end do
                        end do
                        !$omp end do
                else
                        !$omp do
                        do j=1,n
                                do i=1,n
                                        t0 = x(i, j)
                                        t1 = y(i, j)
                                        y(i, j) = t1 * d1(j) + t0 * d2(j)
                                        x(i, j) = t1 * d3(j) + t0 * d1(j)
                                end do
                        end do
                        !$omp end do
                end if
                !$omp end parallel
        end subroutine
        
        subroutine hyperbolic_exp4r(h, n, a, lda, q, r)
                ! real h
                implicit none
                integer, intent(in) :: n, lda
                real(8), intent(in) :: h
                complex(8), intent(in) :: a(lda, n)
                complex(8), intent(out) :: q(lda, n), r(lda, n)
                
                integer :: i
                
                call zgemm('N', 'N', n, n, n, (1d0,0d0), a, lda, a, lda, (0d0,0d0), r, lda)
                q(1:n, :) = ((h**3)/6) * a(1:n, :) + ((h**5)/120) * r(1:n, :)
                do i=1,n
                        q(i, i) = q(i, i) + h
                end do
                r(1:n, :) = (h/2) * a(1:n, :) - ((h**3)/24) * r(1:n, :)
                return
        end subroutine
        
        subroutine hyperbolic_exp6r(h, n, a, q, r)
                implicit none
                integer, intent(in) :: n
                real(8), intent(in) :: h
                complex(8), intent(in) :: a(n, n)
                complex(8), intent(out) :: q(n, n), r(n, n)
                complex(8) :: t(n, n)
                
                real(8), parameter :: cq0 = 1d0/6, cq1 = 1d0/120, cq2 = 1d0/5040, &
                cr0 = 0.5d0, cr1 = -1d0/24, cr2 = 1d0/240
                
                integer :: i
                
                call zgemm('N', 'N', n, n, n, (1d0,0d0), a, n, a, n, (0d0,0d0), r, n)
                call zgemm('N', 'N', n, n, n, (1d0,0d0), r, n, a, n, (0d0,0d0), t, n)
                q = ((h**3)*cq0) * a + ((h**5)*cq1) * r + ((h**7)*cq2) * t
                do i=1,n
                        q(i, i) = q(i, i) + h
                end do
                r = (h*cr0) * a + ((h**3)*cr1) * r + ((h**5)*cr2) * t
                return
        end subroutine
        
        subroutine hyperbolic_exp8cr(h, n, a, q, r)
                implicit none
                integer, intent(in) :: n
                real(8), intent(in) :: h
                complex(8), intent(in) :: a(n, n)
                complex(8), intent(out) :: q(n, n), r(n, n)
                complex :: ca(n, n), t0(n, n), t1(n, n)
                
                integer :: i
                ca = cmplx(a)
                
                q = ((h**3)/6) * a
                do i=1,n
                        q(i, i) = q(i, i) + h
                end do
                r = (h/2) * a
                
                call cgemm('N', 'N', n, n, n, (1.,0.), ca, n, ca, n, (0.,0.), t0, n)
                q = q + ((h**5)/120) * t0
                r = r - ((h**3)/24) * t0
                call cgemm('N', 'N', n, n, n, (1.,0.), ca, n, t0, n, (0.,0.), t1, n)
                q = q + ((h**7)/5040) * t1
                r = r + ((h**5)/240) * t1
                call cgemm('N', 'N', n, n, n, (1.,0.), ca, n, t1, n, (0.,0.), t0, n)
                q = q + ((h**9)/362880) * t0
                r = r - ((h**7)/40320*17) * t0
                
                return
        end subroutine
        
        
        subroutine hyperbolic_exp_adaptiver(h, n, a, q, r, nr0, thre)
                implicit none
                integer, intent(in) :: n
                real(8), intent(in) :: h
                complex(8), intent(in) :: a(n, n)
                complex(8), intent(out) :: q(n, n), r(n, n)
                real(8), intent(in) :: nr0, thre
                
                real(8),parameter,dimension(8) :: cr = [0.5d0, -1d0/24, 1d0/240, -17d0/40320, 31d0/725760, &
                -691d0/159667200d0, 5461d0/12454041600d0, -929569d0/20922789888000d0], &
                cq = [1d0/6, 1d0/120, 1d0/5040, 1d0/362880, 1d0/39916800, &
                1d0/6227020800d0, 1d0/1307674368000d0, 1d0/355687428096000d0]
                
                complex(8) :: t(n, n), t2(n, n)
                integer :: i
                complex(8) :: hh
                real(8) :: e
                hh = h**3
                q = hh * cq(1) * a
                do i=1,n
                        q(i, i) = q(i, i) + h
                end do
                r = h * cr(1) * a
                t = a
                do i=2,7
                        e = abs(((h**2)*hh * cq(i) + hh * abs(cr(i)))) * (nr0**i)
                        if(e.le.thre) exit
                        call zgemm('N', 'N', n, n, n, (1d0,0d0), t, n, a, n, (0d0,0d0), t2, n)
                        t = t2
                        q = q + (h**2) * hh * cq(i) * t
                        r = r + hh * cr(i) * t
                        hh = hh * (h**2)
                end do
                !print *, i, e
                return
        end subroutine
        
        function cond_est(n, A)
                ! use Geshgorin circle
                implicit none
                integer, intent(in) :: n
                complex(8), intent(in) :: A(n, n)
                real(8) :: cond_est
                
                integer :: i, j
                real(8) :: mn, mx, t
                
                mn = abs(A(1,1))
                mx = mn
                !$omp parallel do private(i,j,t) reduction(min:mn) reduction(max:mx)
                do j=1,n
                        t = 0d0
                        do i=1,n
                                if(i.ne.j) then
                                        t = t + sqrt(real(A(i, j))**2 + aimag(A(i,j))**2)
                                end if
                        end do
                        mn = min(mn, (abs(A(j,j)) - t))
                        mx = max(mx, (abs(A(j,j)) + t))
                end do
                !$omp end parallel do
                if(mn < 0d0) then
                        cond_est = 1d100
                else
                        cond_est = mx/mn
                end if
                return
        end function 
        
        function power_est(n, q, p, x, s)
                implicit none
                integer, intent(in) :: n, s
                complex(8), intent(in) :: q(n, n), p(n, n)
                complex(8), intent(inout) :: x(n)
                real(8) :: power_est
                
                complex(8) :: y(n)
                integer :: i
                real(8) :: d0, d, dp
                if(s < 0) then
                        do i=1, n
                                x(i) = ((1d0, 1d0) / sqrt(2d0*n))
                        end do
                        power_est = -1d0
                        return
                end if
                d0 = -1d0
                do i=1,s
                        call power_step(n, q, p, x, y, d, dp)
                        print *, 'P', i, d, dp
                        if(dp < 1d-2 .or. d-d0 < 1d-2*d) exit
                        d0 = d
                end do
                power_est = d
        end function
        
        subroutine power_step(n, q, p, x, z, d, dp)
                implicit none
                integer, intent(in) :: n
                complex(8), intent(in) :: q(n, n), p(n, n)
                complex(8), intent(inout) :: x(n), z(n)
                real(8), intent(inout) :: d, dp
                
                integer :: i
                real(8) :: dn
                complex(8) :: c, cr, y(n+n)
                
                call zgemv('N', n, n, (1d0,0d0), q, n, x, 1, (0d0, 0d0), y, 1)
                call zgemv('N', n, n, (1d0,0d0), p, n, x, 1, (0d0, 0d0), y(n+1), 1)
                call zgemv('C', n, n, (1d0,0d0), q, n, y, 1, (0d0, 0d0), z, 1)
                call zgemv('C', n, n, (1d0,0d0), p, n, y(n+1), 1, (1d0, 0d0), z, 1)
                
                cr = 0d0
                dn = 0d0
                !$omp simd reduction(+:dn,cr) private(c)
                do i=1,n
                        c = z(i)
                        dn = dn + (c%re**2) + (c%im**2)
                        cr = cr + conjg(x(i)) * c
                end do
                !$omp end simd
                d = abs(cr)
                dn = 1d0 / sqrt(dn)
                dp = 0d0
                !$omp simd reduction(+:dp) private(c)
                do i=1, n
                        z(i) = z(i) * dn
                        c = z(i) - x(i)
                        dp = dp + (c%re**2) + (c%im**2)
                        x(i) = z(i)
                end do
                !$omp end simd
                dp = sqrt(dp)
        end subroutine
        
end module
