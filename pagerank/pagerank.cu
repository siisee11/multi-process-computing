/* Each kernel handles the update of one pagerank score. In other
 * words, each kernel handles one row of the update:
 *
 *      pi(t+1) = (1/2) A pi(t) + (1 / (2N))
 *      
 * You may assume that num_nodes <= blockDim.x * 65535
 *
 */
__global__
void device_graph_propagate(const uint* graph_indices
		, const uint* graph_edges
		, const float* graph_nodes_in
		, float* graph_nodes_out
		, const float* inv_edges_per_node
		, int num_nodes) {

	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i < num_nodes){
		float sum = 0.f;

		//for all of its edges
		for(uint j = graph_indices[i]; j < graph_indices[i+1]; j++) {
			sum += graph_nodes_in[ graph_edges[j] ] * inv_edges_per_node[ graph_edges[j] ];
		}

		graph_nodes_out[i] = 0.5f/(float)num_nodes + 0.5f*sum;
	}
}

/* This function executes a specified number of iterations of the
 * pagerank algorithm. The variables are:
 *
 * h_graph_indices, h_graph_edges:
 *     These arrays describe the indices of the neighbors of node i.
 *     Specifically, node i is adjacent to all nodes in the range
 *     h_graph_edges[h_graph_indices[i] ... h_graph_indices[i+1]].
 *
 * h_node_values_input:
 *     An initial guess of pi(0).
 *
 * h_gpu_node_values_output:
 *     Output array for the pagerank vector.
 *
 * h_inv_edges_per_node:
 *     The i'th element in this array is the reciprocal of the
 *     out degree of the i'th node.
 *
 * nr_iterations:
 *     The number of iterations to run the pagerank algorithm for.
 *
 * num_nodes:
 *     The number of nodes in the whole graph (ie N).
 *
 * avg_edges:
 *     The average number of edges in the graph. You are guaranteed
 *     that the whole graph has num_nodes * avg_edges edges.
 *
 */


double device_graph_iterate(const uint* h_graph_indices
		, const uint* h_graph_edges
		, const float* h_node_values_input
		, float* h_gpu_node_values_output
		, const float* h_inv_edges_per_node
		, int nr_iterations
		, int num_nodes
		, int avg_edges) {

	// TODO: allocate GPU memory
	size_t size = num_nodes * sizeof(float);
	float* d_buffer_1;
	cudaMalloc(&d_buffer_1, size);
	float* d_buffer_2;
	cudaMalloc(&d_buffer_2, size);
	uint* d_graph_indices;
	cudaMalloc(&d_graph_indices, (num_nodes+1) * sizeof(int));
	uint* d_graph_edges;
	cudaMalloc(&d_graph_edges, num_nodes * avg_edges * sizeof(int));
	float* d_inv_edges_per_node;
	cudaMalloc(&d_inv_edges_per_node, size);

	// TODO: check for allocation failure

	// TODO: copy data to the GPU
	cudaMemcpy(d_buffer_1, h_node_values_input, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_graph_indices, h_graph_indices, (num_nodes+1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_graph_edges, h_graph_edges, num_nodes * avg_edges * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_inv_edges_per_node, h_inv_edges_per_node , size, cudaMemcpyHostToDevice);

	start_timer(&timer);

	const int block_size = 256;
	const int bnam = num_nodes / block_size + 1;

	
	// TODO: launch your kernels the appropriate number of iterations
	for(int iter = 0; iter < nr_iterations / 2; iter++) {
		device_graph_propagate<<<bnam,block_size>>>(d_graph_indices, d_graph_edges, d_buffer_1, d_buffer_2,
				d_inv_edges_per_node, num_nodes);
		device_graph_propagate<<<bnam,block_size>>>(d_graph_indices, d_graph_edges, d_buffer_2, d_buffer_1,
				d_inv_edges_per_node, num_nodes);
	}
	if(nr_iterations % 2)
		device_graph_propagate<<<bnam,block_size>>>(d_graph_indices, d_graph_edges, d_buffer_1, d_buffer_2,
				d_inv_edges_per_node, num_nodes);

	// This two line below is original code.
	check_launch("gpu graph propagate");
	double gpu_elapsed_time = stop_timer(&timer);

	// TODO: copy final data back to the host for correctness checking
	// handle the odd case and copy memory to the output location
	if(nr_iterations % 2) {
		cudaMemcpy(h_gpu_node_values_output, d_buffer_2, size, cudaMemcpyDeviceToHost);
	} else {
		cudaMemcpy(h_gpu_node_values_output, d_buffer_1, size, cudaMemcpyDeviceToHost);
	}

	// TODO: free the memory you allocated!
	cudaFree(d_buffer_1);
	cudaFree(d_buffer_2);
	cudaFree(d_graph_indices);
	cudaFree(d_graph_edges);
	cudaFree(d_inv_edges_per_node);

	check_launch("gpu cudaFree function");

	return gpu_elapsed_time;
}
