        program overlap_calc
        implicit real*8(a-h,o-z)
        parameter (nmax=25000)
        parameter (nwavetot=10)
        parameter (nmaxtot=nmax*nwavetot)
        dimension x(nmaxtot),weight(nmaxtot),n(nwavetot),n0(nwavetot),
     1  simtime(nwavetot),coord(nmax),coord_sim(nwavetot,nmax),
     1  weight_sim(nwavetot,nmax),tot_weight(nwavetot),reverse_weight(
     1  nwavetot), tot_overlap(nwavetot),ngood(nwavetot)
        open(unit=8,file='testing_d2h.dat',status='old')
        open(unit=9,file='coord_and_weight_dt1_navg560.dat',
     1  status='old')
        open(unit=10,file='overlap_dt1_navg560.dat')
C		Calculates the overlap of two wave functions based on
C		a wave function calculated using DMC. Assumes we are running
C		simulation using importance-sampled DMC.
C		Inputs:
C		coordinates of the walkers and their averaged descendant weights
C		(used 560 for navg)
C		Outputs:
C		The overlap between the wave function calculated using DMC and the
C		wave function used as the trial wave function for DMC. 
        read(9,*) ntot
        read(8,*) nwave
        do k = 1,nwave
            read(8,*) n(k),n0(k),simtime(k)
            do i = 1,n(k)
                read(8,*) coord(i)
            enddo
        enddo
        do i = 1,ntot
            read(9,*) x(i),weight(i)
        enddo
        ip = 0
        do k = 1,nwave
            do i = 1,n(k)
                ip = ip + 1
                coord_sim(k,i) = x(ip)
                weight_sim(k,i) = weight(ip)
            enddo
        enddo
        avg_overlap = 0.d0
        totdev_overlap = 0.d0
        do k = 1,nwave
            tot_weight(k) = 0.d0
            reverse_weight(k) = 0.d0
            ngood(k) = 0
            do i = 1,n(k)
                if (weight_sim(k,i).gt.1e-6) then
                    tot_weight(k) = tot_weight(k) + weight_sim(k,i)
                    reverse_weight(k)=reverse_weight(k)+ 
     1              (1/weight_sim(k,i))
                    ngood(k) = ngood(k) + 1
                endif
            enddo
            tot_overlap(k) = dfloat(ngood(k))/(sqrt(tot_weight(k)*
     1      reverse_weight(k)))
            avg_overlap = avg_overlap + tot_overlap(k)
        enddo
        avg_overlap = avg_overlap/dfloat(nwave)
        do k = 1,nwave 
           totdev_overlap=totdev_overlap+(tot_overlap(k)-avg_overlap)**2
        enddo
        stdev_overlap = sqrt(totdev_overlap/dfloat(nwave))
        write(10,*) avg_overlap,stdev_overlap
        end program     
