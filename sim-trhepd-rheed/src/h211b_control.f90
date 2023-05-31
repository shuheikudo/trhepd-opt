function h211b_limit(u)
        implicit none
        real(8) :: h211b_limit
        real(8), intent(in) :: u
        real(8), parameter :: kappa = 1d0
        h211b_limit = 1d0 + kappa * atan((u-1d0)/kappa)
end function

function h211b_filter(c1, c0, rho, k)
        implicit none
        real(8) :: h211b_filter
        real(8), intent(in) :: c1, c0, rho
        integer, intent(in) :: k
        integer, parameter :: b = 4
        h211b_filter = (c1**(1d0/(b*k))) * (c0**(1d0/(b*k))) * (rho**(-1d0/k))
end function
