!***********************************************************
!       complex matrix  a*inv(b) --> a
!     v.1:84/10  v.4 2014/2/26  T.Hanada
!***********************************************************
subroutine gcmi3(n1,n,a,b)
      implicit none
      integer,intent(in) :: n1,n
      complex(8), intent(inout) :: a(n1, n), b(n1, n)
      complex(8) :: t(n, n)
      integer :: p(n), info
      call zgetrf(n, n, b, n1, p, info)
      t = transpose(a(1:n, 1:n))
      call zgetrs('T', n, n, b, n1, p, t, n, info)
      a(1:n, 1:n) = transpose(t)
end subroutine