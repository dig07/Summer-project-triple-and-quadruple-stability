
subroutine simulate(i,a_ratio1,a_ratio2,qoutlist,size_qoutlist_0,size_qoutlist_1,samples_eccentricities, &
	size_samples_eccentricities,inclinations, size_inclinations,stable_aratio_counter,&
	primary_eccentricities_counter,second_eccentricities_counter,qout_counter)
	
	integer*8 :: stable_ain_counter,stable_qout_counter,stable_qout_counter_primary,stable_qout_counter_secondary
	integer :: i,size_aout,size_qoutlist,size_samples_eccentricities,size_inclinations
	integer:: size_qoutlist_0,size_qoutlist_1
	real, dimension(:), allocatable :: stable_ecc,stable_ecc_2

	real ::a_ratio1,a_ratio2
	real, dimension(size_qoutlist_0,size_qoutlist_1) :: qoutlist
        real, dimension(size_samples_eccentricities) :: samples_eccentricities
	real, dimension(size_inclinations) :: inclinations
	double precision :: inclination_limit,inclination_limit_2

	integer*8, dimension(size_qoutlist_0,size_qoutlist_1)::qout_counter
	integer*8, dimension(size_samples_eccentricities)::primary_eccentricities_counter,second_eccentricities_counter
	real:: Pi
	Pi = 3.1415926535
	!f2py intent(in) :: i,a_ratio1,a_ratio2,qoutlist,samples_eccentricities,inclinations
	!f2py intent(hide), depend(qoutlist) :: size_qoutlist_0= shape(qoutlist,0)
	!f2py intent(hide), depend(qoutlist) :: size_qoutlist_1= shape(qoutlist,1)
	!f2py intent(hide), depend(samples_eccentricities) :: size_samples_eccentricities= shape(samples_eccentricities,0)
	!f2py intent(hide), depend(inclinations) :: size_inclinations= shape(inclinations,0)
	!f2py intent(out) :: stable_aratio_counter,primary_eccentricities_counter,second_eccentricities_counter,qout_counter

	stable_aratio_counter = 0
	qout_counter = 0
	eccentricities_counter = 0
	primary_eccentricities_counter = 0 
	secondary_eccentricities_counter = 0 

	do k=1,size(qoutlist(1,:))
		qout = qoutlist(1,k)
		stable_qout_counter_primary = 0
		do m=1,size(samples_eccentricities)
			eccentricity = samples_eccentricities(m)
			primary_eccentricity_counter = 0
			inclination_limit = (Pi/0.3)*(1-(a_ratio1)/((( (1+qout)*(1+eccentricity)/(sqrt(1-& 
			      eccentricity)))**(2./5.))*(2.8/(1-eccentricity))))
			stable_ecc = pack(inclinations,inclinations>inclination_limit)
			stable_1 = int(size(stable_ecc))
			if (stable_1.NE.0) then
				do z =2,size(qoutlist(:,k))
					qout_secondary = qoutlist(z,k)
					stable_qout_counter_secondary = 0
					do r=1,size(samples_eccentricities)
						secondary_eccentricity = samples_eccentricities(r)
						inclination_limit_2 = (Pi/0.3)*(1-(a_ratio2)/((( (1+qout_secondary)*(1+secondary_eccentricity)/(sqrt(1-& 
							secondary_eccentricity)))**(2./5.))*(2.8/(1-secondary_eccentricity))))						
						stable_ecc_2 = pack(inclinations,inclinations>inclination_limit_2)
						stable_2 = int(size(stable_ecc_2))
						stable_qout_counter_primary = stable_qout_counter_primary + stable_2*stable_1
						stable_qout_counter_secondary = stable_qout_counter_secondary +stable_2*stable_1
						secondary_eccentricities_counter(r) = secondary_eccentricities_counter(r) + stable_2*stable_1
						primary_eccentricity_counter = primary_eccentricity_counter + stable_2*stable_1
						stable_aratio_counter = stable_aratio_counter +stable_2*stable_1
					end do
					qout_counter(z,k)=qout_counter(z,k)+stable_qout_counter_secondary
				end do
			end if
			primary_eccentricities_counter(m) = primary_eccentricities_counter(m)+primary_eccentricity_counter
		end do
		qout_counter(1,k)=qout_counter(1,k)+stable_qout_counter_primary
	end do
end subroutine simulate 
