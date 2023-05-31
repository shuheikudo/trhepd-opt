!*******************************************************************
!   rheed multi slice method
!   subroutine trmat gcmi2
!   v.1:84/10   v.2:86/11   v.3:90/4   v.4:2014/4   v.4b:2017/2
!    T.Hanada
!*******************************************************************
subroutine srfref(nv,nb,ns,iv,dz,ngr,nbg,iord,ih,ik &
                ,wn,ghx,ghy,gky,azi,daz,naz,gi,dg,ng,idiag,iprn,int_method)
        use mod_scpot3, only: comp_scpot
        use mod_prkn, only : hp_int4, hp_int_init4, hp_int6, hp_int_init6,hp_int6_1d
        use mod_rk, only : rk4_int, rk4_int_init, rk4_intgs
        use mod_magnus, only : m2_int_init, m4_int_init, &
                lf_int_spa, m2_int_spa, m4_int_spa, m6_int_spa

        implicit none
        real(8) :: dz,wn,ghx,ghy,gky,azi,daz,gi,dg
        integer :: iv(nb,nb),ih(nb),ik(nb),iord(nb),nbg(ngr)
        integer :: nv,nb,ns,ngr,naz,ng,idiag,iprn, int_method

        complex(8) :: t1(nb,nb,max(ng,naz)),gma(nb)
        complex(8),allocatable :: f(:, :)
        real(8) :: f2(nb, max(ng, naz))
        integer :: nb2,nrep1,nrep2,irep1,irep2,nb0,ib1,ib2,i,j,l
        real(8) :: az,caz,saz,ga,cga,sga,angle,wnx,wny,wgx,wgy,wnsga2,s


        complex(8), allocatable, dimension(:,:) :: v, vi
        select case (int_method)
        case (1,4)
                call m2_int_init(dz, ns, nv, v, vi)
        case (2:3)
                call m4_int_init(dz, ns, nv, v, vi)
        case (5)
                call hp_int_init4(dz, ns, nv, v, vi)
        case (6)
                call hp_int_init6(dz, ns, nv, v, vi)
        case (7,8)
                call rk4_int_init(dz, ns, nv, v, vi)
        case default
                call abort
        end select

        if (ng > naz) then
                nrep1=naz
                nrep2=ng
        else
                nrep1=ng
                nrep2=naz
        endif
        nb2=nb+nb
        do i=1,nb
                if (ih(i) == 0 .and. ik(i) == 0) nb0=i
        end do
        !----------azimuth,glancing angle scan----------
        do irep1=1,nrep1
                if (ng > naz) then
                        az=azi+(irep1-1)*daz
                        caz=cos(az)
                        saz=sin(az)
                else
                        ga=gi+(irep1-1)*dg
                        cga=cos(ga)
                        sga=sin(ga)
                endif
                do irep2=1,nrep2
                        ib1=1
                        do l=1,ngr
                                ib2=ib1+nbg(l)-1
                                read (1) ((t1(i,j,irep2),i=ib1,ib2),j=ib1,ib2)
                                ib1=ib2+1
                        end do
                end do
                !$omp parallel do &
                !$omp firstprivate(ga,angle,cga,sga,az,caz,saz)&
                !$omp private(wnsga2,wnx,wny,wgx,wgy,s,gma,f,i,&
                !$omp   ib1,l,ib2,j,irep2) &
                !$omp if(nb<100 .and. nrep2>20) &
                !$omp schedule(dynamic)
                do irep2=1,nrep2
                        allocate(f(nb,nb))
                        if (ng > naz) then
                                ga=gi+(irep2-1)*dg
                                cga=cos(ga)
                                sga=sin(ga)
                        else
                                az=azi+(irep2-1)*daz
                                caz=cos(az)
                                saz=sin(az)
                        endif
                        wnsga2=wn*sga*sga
                        wnx=wn*cga*caz
                        wny=wn*cga*saz
                        !----------z component of scattering vector----------
                        do i=1,nb
                                wgx=ghx*ih(i)+wnx
                                wgy=ghy*ih(i)+gky*ik(i)+wny
                                s=wn*wn-wgx*wgx-wgy*wgy
                                if (s >= 0d0) then
                                        gma(i)=dcmplx(sqrt(s))
                                else
                                        gma(i)=dcmplx(0d0,sqrt(-s))
                                endif
                        end do
                        !----------diffraction from bulk layer----------
                        f(:,:)=(0d0,0d0)
                        ib1=1
                        do l=1,ngr
                                ib2=ib1+nbg(l)-1
                                do  j=ib1,ib2
                                        do  i=ib1,ib2
                                                f(iord(i),iord(j))=t1(i,j,irep2)
                                        end do
                                end do
                                ib1=ib2+1
                        end do

                        ! integration
                        select case (int_method)
                        case (1)
                                call m2_int_spa(dz, ns, nb, nv, f, gma, iv, v, vi)
                        case (2)
                                call m4_int_spa(dz, ns, nb, nv, f, gma, iv, v, vi)
                        case (3)
                                call m6_int_spa(dz, ns, nb, nv, f, gma, iv, v, vi)
                        case (4)
                                call lf_int_spa(dz, ns, nb, nv, f, gma, iv, v, vi)
                        case (5)
                                call hp_int4(dz, ns, nb, nv, f, gma, iv, v, vi)
                        case (6)
                                if(nb .eq. 1) then
                                        call hp_int6_1d(dz, ns, nb, nv, f, gma, iv, v, vi)
                                else
                                        call hp_int6(dz, ns, nb, nv, f, gma, iv, v, vi)
                                end if
                        case (7)
                                call rk4_int(dz, ns, nb, nv, f, gma, iv, v, vi)
                        case (8)
                                call rk4_intgs(dz, ns, nb, nv, f, gma, iv, v, vi)
                        end select


                        !----------diffracion intensity---------
                        do i=1,nb
                                s=dble(gma(i))
                                if (s > 1d-10) then
                                        ! wnsga2=wn*sga*sga=K*sin^2(theta_g)=gma(nb0)*sin(theta_g)
                                        f2(i,irep2)=dble(f(i,nb0)*dconjg(f(i,nb0)))*wnsga2/s
                                else
                                        f2(i,irep2)=0d0
                                endif
                        end do
                        deallocate(f)
                end do ! irep2=1,nrep2
                do irep2=1,nrep2
                        if (ng > naz) then
                                ga=gi+(irep2-1)*dg
                                angle=ga*180d0/3.141592654d0
                        else
                                az=azi+(irep2-1)*daz
                                angle=az*180d0/3.141592654d0
                        endif
                        write (3,'(E12.4,200(",",E23.15))') angle,(f2(i,irep2),i=1,nb)
                        ! print ('(E12.4,200(",",E12.4))'), angle,(abs(f(i,nb0)),i=1,nb)
                end do
                write (3,*)
        end do ! irep1=1,nrep1

        deallocate(vi)
        deallocate(v)
        return
end subroutine


! use leapforg integration
! about ~10x faster, with O(10^{-6}) abs-error
! compared with surf.f90
subroutine lf_int(dz, ns, nb, nv, f, g, iv, v, vi)
        use mod_matcomp
        implicit none
        integer, intent(in) :: ns, nb, nv, iv(nb,nb)
        real(8), intent(in) :: dz
        complex(8), intent(in) :: g(nb), v(nv, ns), vi(nv, ns)
        complex(8), intent(inout) :: f(nb,nb)

        complex(8) :: a, v0(nv), v1(nv), t1(nb, nb+nb)
        integer :: i, j, l
        complex(8) :: xz, xz2

        xz = (0d0,-1d0)*dz
        xz2 = (0d0,-0.5d0)*dz

        ! transform variables to reconstruct the simplest form
        ! \"{y} = f(y).
        ! F <- (I-F) / (I+F) * G
        do j=1,nb
                do i=1,nb
                        if( i .eq. j) then
                                t1(i, j) = 1d0 - f(i, j)
                                t1(i, nb+j) = 1d0 + f(i, j)
                        else
                                t1(i, j) = -f(i, j)
                                t1(i, nb+j) = f(i, j)
                        end if
                end do
        end do
        call matAdivB(nb, t1, nb)
        do j=1,nb
                do i=1,nb
                        t1(i, nb+j) = t1(i, j) * g(j)
                end do
        end do
        ! in the rest of the process, the matrix is located at t1(i, nb+j)


        ! the leapfrog integration
        ! exp(h*[0 U(z); 1 0]) \approx
        ! [1 0; h/2 1][1 \int_z^{z+h}U(z)dz; 0 1][1 0; h/2 1]
        ! The first q integration in inverse form:
        ! F <- (-0.5i * dz * F  + I) \ F
        do j=1,nb
                do i=1,nb
                        if(i .eq. j) then
                                t1(i, j) = xz2*t1(i, nb+j) + 1d0
                        else
                                t1(i, j) = xz2*t1(i, nb+j)
                        end if
                end do
        end do
        call matAbdivB(nb, t1, nb)

        do l=1,ns
                ! p integration
                ! F <- F + -1i*dz*(V+Vi+G*G)
                v0 = v(:,l) + vi(:,l)
                v1 = dconjg(v(:,l)-vi(:,l))

                do j=1,nb
                        do i=1,j-1
                                t1(i, nb+j) = t1(i, nb+j) + xz*v1(iv(j,i))
                        end do
                        t1(j, nb+j) = t1(j, nb+j) + xz*(v0(1) + g(j)*g(j)) 
                        do i=j+1,nb
                                t1(i, nb+j) = t1(i, nb+j) + xz*v0(iv(i,j))
                        end do
                end do

                if(l .eq. ns) exit ! skip last q int.

                ! q integration
                ! F <- (-1i*dz*F+I) \ F
                do j=1,nb
                        do i=1,nb
                                if(i .eq. j) then
                                        t1(i, j) = 1d0 + xz*t1(i, nb+j)
                                else
                                        t1(i, j) = xz*t1(i, nb+j)
                                end if
                        end do
                end do
                call matAbdivB(nb, t1, nb)
        end do ! l=1,ns
        ! p integration of the last step
        do j=1,nb
                do i=1,nb
                        if(i .eq. j) then
                                t1(i, j) = 1d0 + xz2*t1(i, nb+j)
                        else
                                t1(i, j) = xz2*t1(i, nb+j)
                        end if
                end do
        end do
        call matAbdivB(nb, t1, nb)

        ! back-transform of the variables
        do j=1,nb
                do i=1,nb
                        a = t1(i, nb+j) * (1d0/g(j))
                        if( i .eq. j) then
                                t1(i, j) = 1d0 - a
                                t1(i, nb+j) = 1d0 + a
                        else
                                t1(i, j) = -a
                                t1(i, nb+j) = a
                        end if
                end do
        end do
        call matAdivB(nb, t1, nb)

        f(1:nb, 1:nb) = t1(1:nb, 1:nb)
end subroutine
