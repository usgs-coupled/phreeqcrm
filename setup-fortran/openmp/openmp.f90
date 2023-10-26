program main
use omp_lib

integer threads
integer partial_sum, sum

threads = 0
!$OMP parallel reduction(+ : threads)
threads = threads + 1
!$OMP END PARALLEL

sum = 0
!$OMP DO
do i = 0, 1000
    partial_sum = partial_sum + i
end do
!$OMP END DO

!$OMP CRITICAL
sum = sum + partial_sum
!$OMP END CRITICAL

print *, 'using ', threads, 'threads, sum(0:1000) = ', sum
end program
