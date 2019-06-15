#pragma omp parallel default(private) shared(npoints) reduction(+:sum) num_thread(8)
{
	num_threads = omp_get_num_threads();
	sample_points_per_thread = npoints / num_threads;
	sum = 0;
	int i = 0;
	for (i = 0 ; i < sample_points_per_thread ; i++){
		rand_no_x = (double)rand_r(&seed)/RAND_MAX;
		rand_no_y = (double)rand_r(&seed)/RAND_MAX;
		if ( x * x + y * y < 1.0)
			sum ++;
	}
}

