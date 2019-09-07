
subroutine simulate(i,a_ratio,qoutlist,size_qoutlist,samples_eccentricities,size_samples_eccentricities, &
	inclinations, size_inclinations,stable_aratio_counter,eccentricities_counter,qout_counter)
	integer*8 :: stable_ain_counter,stable_qout_counter
	integer :: i,size_aout,size_qoutlist,size_samples_eccentricities,size_inclinations
	real, dimension(:), allocatable :: stable_ecc
	real ::a_ratio
    	real, dimension(size_qoutlist) :: qoutlist
    	real, dimension(size_samples_eccentricities) :: samples_eccentricities
	real, dimension(size_inclinations) :: inclinations
	double precision :: inclination_limit

	integer*8, dimension(size_qoutlist)::qout_counter
	integer*8, dimension(size_samples_eccentricities)::eccentricities_counter
	real:: Pi
	Pi = 3.1415926535
	!f2py intent(in) :: i,a_ratio,qoutlist,samples_eccentricities,inclinations
	!f2py intent(hide), depend(qoutlist) :: size_qoutlist= shape(qoutlist,0)
	!f2py intent(hide), depend(samples_eccentricities) :: size_samples_eccentricities= shape(samples_eccentricities,0)
	!f2py intent(hide), depend(inclinations) :: size_inclinations= shape(inclinations,0)
	!f2py intent(out) :: stable_aratio_counter,eccentricities_counter,qout_counter

	stable_aratio_counter = 0
	qout_counter = 0
	eccentricities_counter = 0

	do k=1,size(qoutlist)
		qout = qoutlist(k)
		stable_qout_counter = 0
		do m=1,size(samples_eccentricities)
			eccentricity = samples_eccentricities(m)
			inclination_limit = (Pi/0.3)*(1-(a_ratio)/((( (1+qout)*(1+eccentricity)/(sqrt(1-& 
			      eccentricity)))**(2./5.))*(2.8/(1-eccentricity))))
			stable_ecc = pack(inclinations,inclinations>inclination_limit)
			stable_qout_counter = stable_qout_counter + size(stable_ecc)
			stable_aratio_counter = stable_aratio_counter +size(stable_ecc)
			eccentricities_counter(m)= eccentricities_counter(m)+size(stable_ecc)
		end do
		qout_counter(k)=qout_counter(k)+stable_qout_counter
	end do
end subroutine simulate 
