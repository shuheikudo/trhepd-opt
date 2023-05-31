!*******************************************************************
!   rheed multi slice method
!   subroutine trmat gcmi2
!   v.1:84/10   v.2:86/11   v.3:90/4   v.4:2014/4   v.4b:2017/2
!    T.Hanada
!*******************************************************************
       subroutine srfref(nv,nb,ns,v,vi,iv,dz,ngr,nbg,iord,ih,ik &
                  ,wn,ghx,ghy,gky,azi,daz,naz,gi,dg,ng,idiag,iprn)
        implicit none
        complex(8),intent(in) :: v(nv,ns),vi(nv,ns)
        real(8),intent(in) :: dz,wn,ghx,ghy,gky,azi,daz,gi,dg
        integer,intent(in) :: iv(nb,nb),ih(nb),ik(nb),iord(nb),nbg(ngr)
        integer,intent(in) :: nv,nb,ns,ngr,naz,ng,idiag,iprn

        complex(8),allocatable :: t2(:,:), t3(:,:), f(:,:)
        complex(8) :: gma(nb)
        complex(8) :: t1(nb+nb, nb+nb, max(naz,ng))
        real(8) :: f2(nb, max(naz,ng))
        integer :: nb2,nrep1,nrep2,irep1,irep2,nb0,ib1,ib2,i,j,j2,l,k
        real(8) :: az,caz,saz,ga,cga,sga,angle,wnx,wny,wgx,wgy,wnsga2,s

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
          do irep2=1, nrep2
            ib1=1
            do l=1,ngr
              ib2=ib1+nbg(l)-1
              read (1) ((t1(i,j,irep2),i=ib1,ib2),j=ib1,ib2)
              ib1=ib2+1
            end do
          end do

          !$omp parallel do &
          !$omp firstprivate(ga,angle,cga,sga,az,caz,saz) &
          !$omp private(wnsga2,wnx,wny,wgx,wgy,s,gma,i,f, &
          !$omp   ib1,l,ib2,j,irep2,t2,t3) &
          !$omp if(nb<100 .and. nrep2>20) &
          !$omp schedule(static)
          do irep2=1,nrep2
            allocate(t2(nb+nb,nb+nb))
            allocate(t3(nb+nb,nb+nb))
            allocate(f(nb,nb))
            if (ng > naz) then
              ga=gi+(irep2-1)*dg
              angle=ga*180d0/3.141592654d0
              cga=cos(ga)
              sga=sin(ga)
            else
              az=azi+(irep2-1)*daz
              angle=az*180d0/3.141592654d0
              caz=cos(az)
              saz=sin(az)
            endif
            wnsga2=wn*sga*sga
            wnx=wn*cga*caz
            wny=wn*cga*saz
            if (iprn == 1) write (*,*) ' irep1,irep2=',irep1,irep2
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
            !----------transfer matrix of 1 slice (t1)-----------
            do l=1,ns
              call trmat(nv,nb,nb,gma,v(1,l),vi(1,l),iv &
                        ,dz,t1(1,1,irep2),t2,t3,idiag)
              !------------diffraction from surface layer----------
              t2(1:nb,1:nb) = t1(1:nb,(nb+1):(2*nb),irep2)
              t3(1:nb,1:nb) = t1((nb+1):(2*nb),(nb+1):(2*nb),irep2)
              call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1(1,1,irep2), 2*nb, f, nb, (1d0,0d0), t2, 2*nb)
              call zgemm('N', 'N', nb, nb, nb, (1d0,0d0), t1(nb+1,1,irep2), 2*nb, f, nb, (1d0,0d0), t3, 2*nb)
              call gcmi3(nb2,nb,t2,t3)
              f(1:nb,1:nb)=t2(1:nb,1:nb)
            end do ! l=1,ns
            !----------diffracion intensity---------
            do i=1,nb
              s=dble(gma(i))
              if (s > 1d-10) then
                ! wnsga2=wn*sga*sga=K*sin^2(theta_g)=gma(nb0)*sin(theta_g)
                f2(i, irep2)=dble(f(i,nb0)*dconjg(f(i,nb0)))*wnsga2/s
              else
                f2(i, irep2)=0d0
              endif
            end do
            deallocate(f)
            deallocate(t3)
            deallocate(t2)
          end do ! irep2=1,nrep2
          do irep2=1,nrep2
            if (ng > naz) then
              ga=gi+(irep2-1)*dg
              angle=ga*180d0/3.141592654d0
            else
              az=azi+(irep2-1)*daz
              angle=az*180d0/3.141592654d0
            endif
            write (3,'(E12.4,200(",",E23.15))') angle,(f2(i, irep2),i=1,nb)
          end do
          write (3,*)
        end do ! irep1=1,nrep1

        return
        end
