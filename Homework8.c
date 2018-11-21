#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>




// SIMULATION PARAMETERS
#define SIM_TIME 5000 // Simulation time
#define dt 0.00001 // Time Increment
#define N 100 // Number of particles
#define INIT_SPEED 10.0 // Maximum inital speed of particles
#define FILENAME "ParticleVelocity.dat"


// Other params to fine tune simulation
#define epsilon 1.0
#define sigma 1.0
// Dimensional quatities

#define L_d 10 // box sze L_d x L_d
#define beta 1.0 // 1/boltzmann*T 
#define particle_mass 1.0
// Boundary walls threshold set tp 0 for no padding
#define threshold 0.0


/*

Particle simulated using C struct with following fields.

(Double)
pos_x : x coordinate of position
pos_y : y coordinate of position

vel_x : x component of velocity
vel_y : y component of velocity


acc_x : x component of acceleration
acc_y : y component of acceleration

acc_x_next_step : x component of acceleration at next timestep
acc_y_next_step : y component of acceleration at next timestep

(int)
id : particle identification used to ensure a particle is not measured against itself

*/
struct Particle {
	// In two dimensions
	int id; // particle identification
	double pos_x;
	double pos_y;
	double vel_x;
	double vel_y;
	double acc_x;
	double acc_y;
	// acceleration for next time step
	double acc_x_next_step;
	double acc_y_next_step;
	// Mass assumed constant over all particles
};


/*

Particle Functions

*/

// Test if particle within bounds of box
int particle_test(struct Particle p){
	if( (p.pos_x > 0 && p.pos_x < L_d) && (p.pos_y > 0 && p.pos_y < L_d) ){
		return 0; // within box
	} else {
		return 1;
	}
}

// Distance between two points in cartesian space
double particle_separation_distance(struct Particle p1, struct Particle p2){
	return sqrt(pow((p1.pos_x - p2.pos_x),2.0) + pow((p1.pos_y - p2.pos_y),2.0));
}

// Calculate magnitude of vector given x and y components
double magnitude(double x, double y){
	return sqrt(pow(x,2) + pow(y,2));
}

// Given intermolecular seperation between two particles return force between them
double force(double r){
	return -4.0*epsilon*( (-12.0*pow(sigma,12)*pow(r,-11)) + (6.0*pow(sigma,6) * pow(r,-5)) );
}

// Get force in x direction due to interaction of two particles
double force_x(struct Particle p1, struct Particle p2){
	return force(particle_separation_distance(p1,p2)) * ( (p1.pos_x - p2.pos_x) / sqrt(pow(p1.pos_x - p2.pos_x,2.0) + pow(p1.pos_y - p2.pos_y, 2.0)) );
}

// Get force in x direction due to interaction of two particles
double force_y(struct Particle p1, struct Particle p2){
	return force(particle_separation_distance(p1,p2)) * ( (p1.pos_y - p2.pos_y) / sqrt(pow(p1.pos_x - p2.pos_x,2.0) + pow(p1.pos_y - p2.pos_y, 2.0)) );
}

// Reads in a number and returns the opposite sign of the number
double flip_sign(double x){
	if (x > 0){
		return -1.0;
	} else{
		return 1.0;
	}
}

// Update velocity of particle at boundary as given by [van Beijeren, 2014, http://arxiv.org/abs/1411.2983]
double update_velocity(double old_vel){
	return flip_sign(old_vel) * sqrt( (-2.0/(beta*particle_mass)) *  log( 1.0 - exp((-0.5 * beta * particle_mass) * pow(old_vel,2.0))) );
}

// returns the number of particles that left the box at a given timestep
int misbehaving_particle_counter(struct Particle p_array[N]){
	int count = 0;
	for(int i = 0; i < N; i++){
		if(particle_test(p_array[i])) {
			count++;
			printf("Particle ID %2.1d | pos: (%f,%f) | vel: (%f,%f) | acc: (%f,%f)\n",p_array[i].id,p_array[i].pos_x,p_array[i].pos_y,p_array[i].vel_x,p_array[i].vel_y,p_array[i].acc_x,p_array[i].acc_y);
		}
	}
	return count;
}


/*

MISC Functions 

*/


// return random number within dimensions of box
double rand_pos_component() {
    return ((double)L_d)*(double)rand() / (double)RAND_MAX ;
}

// return random number in range [0,INIT_SPEED]
double rand_val_component() {
	return (double)(INIT_SPEED)*((double)rand() / (double)RAND_MAX) ;
}


// check if value v in array up to z
int checkArray(double arr[N], double v, int z){
	for(int i = 0; i < z; i++){
		if(arr[i] >= v){
			return 1; // found in array
		}
	}
	return 0; // not found in array
}

// Print out data on all particles
void particlePrinter(struct Particle p_array[N]){
	for(int i = 0; i < N; i++){
		printf("Particle ID %2.1d | pos: (%f,%f) | vel: (%f,%f) | acc: (%f,%f)\n",p_array[i].id,p_array[i].pos_x,p_array[i].pos_y,p_array[i].vel_x,p_array[i].vel_y,p_array[i].acc_x,p_array[i].acc_y);
	}
}


/*Write out data to a .dat file for plotting*/
int WriteToFile(const char *Filename, int size, double *J){
    int i;
    FILE *fptr;
    char there_was_error = 0;

    /* File Handling */
    fptr = fopen(Filename, "w+");
    if(fptr == NULL) //if file does not exist, create it
    {

        fptr = fopen(Filename, "a+");
        if (fptr == NULL)
            there_was_error = 1;
    }

    if (there_was_error)
    {
        printf("Disc full or no permission\n");
        return 0;
    }

    /* File Writing */

  /*Two column writing line by line per iteration*/
    for(i=0; i< size; i++) {
        fprintf(fptr,"%d    %*.*f\n",i,5,9,J[i]);
    }
	/* Close file after write out*/
    fclose(fptr);
    return 1;
}


/*

Particle Energy Functions

*/

// Potential energy
double potential(struct Particle p, struct Particle p_array[N]){
	// Loop over all particles
	// r seperation distance
	double r;
	double energy = 0;
	for(int i = 0; i < N; i++){
		// Compute force due to all other particles
		if(p.id != p_array[i].id){
			r = particle_separation_distance(p,p_array[i]);

			energy += 4.0*epsilon*( pow((sigma/r),12.0) - pow((sigma/r),6.0) );
		}
	}

	return (energy);
}

// Kinetic energy ( 0.5mv^2 )
double kinetic(struct Particle p){
	return 0.5 * pow(magnitude(p.vel_x,p.vel_y),2.0);
}

// Calculate total energy of single particle
double total_energy(struct Particle p, struct Particle p_array[N]){
	return potential(p,p_array) + kinetic(p);
}

/*

Statistics Functions

- Used to compute different statistics for simulation


*/



/*

We only need distances from the * elements in the following incidence matrix of particles. (x is diagonal).

|		x
|     x *
|   x * *
| x * * *
 __________


*/
double Avg_Seperation(struct Particle p_array[N]){
	int count = 0.0;
	double sep = 0.0;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < i ; j++) {
			count+=1.0;
			sep+=particle_separation_distance(p_array[i],p_array[j]);
		}
	}
	return sep/count;
} // Avg_Seperation per timestep

double Internal_energy(struct Particle p_array[N]) {
	double agg_energy = 0.0;
	for(int i = 0; i < N; i++){
		agg_energy+=total_energy(p_array[i],p_array);
	}
	return agg_energy;
} // Internal_energy per timestep

double Avg_velocity(struct Particle p_array[N]){
	double agg_speed = 0.0;
	for(int i = 0; i < N; i++) {
		agg_speed += magnitude(p_array[i].vel_x,p_array[i].vel_y);
	}
	return agg_speed/(double)N;
} // Avg_velocity per timestep




/*

Particle Boundary Condition Functions

*/

// If the flip_on_boundary causes a particle to be placed on top of another particle this function will correct that particle
int particle_correction(struct Particle p_array[N], int id){
	double sep = 0.00000001;
	int ret = 0;
	for(int i = 0; i < N; i++){
		if( i != id){
			while(particle_separation_distance(p_array[i], p_array[id]) < sep){
				//printf("PARTICLE CORRECTION REQUIRED (%f %f)\n ", p_array[i].pos_x, p_array[i].pos_y);
				if(p_array[id].pos_x < 5.0){
					p_array[id].pos_x += sep*10.0;
				} else {
					p_array[id].pos_x -= sep*10.0;
				}

				if(p_array[id].pos_y < 5.0){
					p_array[id].pos_y += sep*10.0;
				} else {
					p_array[id].pos_y -= sep*10.0;
				}


				ret = 1;	
			}
		}
	}
	return ret;
}


void flip_on_boundary(struct Particle *p) {

	/*
	Example

	T0 					T1

	* |   		->   	| *
	  |					|

	*/

	if( p->pos_x <= 0.0) {
		if (p->pos_x > -1.0*L_d){
		p->pos_x = -1.0 * p->pos_x;
		} else{
			p->pos_x = 1.0;
		}

	}

	if( p->pos_x >= L_d) {
		if(p->pos_x < (2.0*L_d)) {
		p->pos_x = L_d - (p->pos_x - L_d);
		} else { // Out of bounds by 2*L_d
			p->pos_x = L_d - 1.0;
		}
	}

	if( p->pos_y <= 0.0) {
		if (p->pos_y > -1.0*L_d){
		p->pos_y = -1.0 * p->pos_y;
		} else{
			p->pos_y = 1.0;
		}
	}
	if( p->pos_y >= L_d) {
		if(p->pos_y < (2.0 * L_d)){
		p->pos_y = L_d - (p->pos_y - L_d);
		} else { // Out of bounds by 2*L_d
			p->pos_y = L_d - 1.0;
		}
	}

}


void enforce_boundary(struct Particle p_array[N]){
	// Loop through all particles and check if they are near boundary.
	// If near boundary then particle needs to rotate by 180 degrees
	for(int i = 0; i < N; i++){
		// If particle hits left or right side of box where position of x coordinate near 0 or L
		if( (p_array[i].pos_x <= threshold) || (p_array[i].pos_x >= (L_d - threshold)) ){
			// Keep y component of velocity and update x component
			p_array[i].vel_x = update_velocity(p_array[i].vel_x);

		} if( (p_array[i].pos_y <= threshold) || (p_array[i].pos_y >= (L_d - threshold)) ){
			// Keep x component of velocity and update y component
			p_array[i].vel_y = update_velocity(p_array[i].vel_y);

		}
	}
}

// Enforce boundary cond on  a single particle PASS BY POINTER
void enforce_boundary_individual_particle(struct Particle p){
	// Loop through all particles and check if they are near boundary.
	// If near boundary then particle needs to rotate by 180 degrees

	// If particle hits left or right side of box where position of x coordinate near 0 or L
	if( (p.pos_x <= threshold) || (p.pos_x >= (L_d - threshold)) ){
		// Keep y component of velocity and update x component
		p.vel_x = update_velocity(p.vel_x);

	} if( (p.pos_y <= threshold) || (p.pos_y >= (L_d - threshold)) ){
		// Keep x component of velocity and update y component
		p.vel_y = update_velocity(p.vel_y);

	}

}


/*

Particle Motion Functions
- Velocity Verlet Algorithm 
- Update Acceleraion due to forces of other particles

*/

// Update Acceleration of ALL particles
void UpdateAcceleration(struct Particle p_array[N]){
	double f_x, f_y;
	for(int i =0; i < N; i++){ 	// loop though each particle
		f_x = 0.0;
		f_y = 0.0;
		for(int j = 0; j < N; j++){ // Calculate force of all other particles on ith particle
			if(p_array[i].id != p_array[j].id){
				f_x+=force_x(p_array[i],p_array[j]);
				f_y+=force_y(p_array[i],p_array[j]);
			}
		}
		// Update particle

		p_array[i].acc_x = p_array[i].acc_x_next_step;
		p_array[i].acc_y = p_array[i].acc_y_next_step;
		p_array[i].acc_x_next_step = f_x;
		p_array[i].acc_y_next_step = f_y;
	}
}
// Velocity Verlet Algorithm position update step ** See documentation **
// Will update positions of all particles in accordance with algorithm and store results in each particles struct
void UpdatePosition(struct Particle p_array[N]){
	for(int i = 0; i < N; i++){
		// We update position value of particle
		double x_old = p_array[i].pos_x;
		double y_old = p_array[i].pos_y;

		p_array[i].pos_x = p_array[i].pos_x + p_array[i].vel_x*dt + 0.5*p_array[i].acc_x*pow(dt,2.0);
		p_array[i].pos_y = p_array[i].pos_y + p_array[i].vel_y*dt + 0.5*p_array[i].acc_y*pow(dt,2.0);

		// If particle has gone out bounds as a result of the update position

		if(particle_test(p_array[i])) {

			//printf("HIT CONDITION %d \n",i );
			// Update velocity
			enforce_boundary_individual_particle(p_array[i]);

			// reflect particle position within box
			flip_on_boundary(&p_array[i]);

			while(particle_correction(p_array,i)); // Ensure particle isnt on top and correct if it is
		}
	}
}

// Velocity Verlet Algorithm velocity update step ** See documentation **
// Will update velocity of all particles in accordance with algorithm and store results in each particles struct
void UpdateVelocity(struct Particle p_array[N]){
	for(int i = 0; i < N; i++){
		p_array[i].vel_x = p_array[i].vel_x + 0.5*dt*(p_array[i].acc_x + p_array[i].acc_x_next_step);
		p_array[i].vel_y = p_array[i].vel_y + 0.5*dt*(p_array[i].acc_y + p_array[i].acc_y_next_step);
	}
}







/*
Particle Intialization methods

*/

/*
The following method will randomly assign locations to particles.
** MAY CAUSE PARTICLES TO GENERATE WITHIN CLOSE PROXIMITY **
*/
void initialize_particles(struct Particle p_array[N]){
	double init_x_pos[N];
	double init_y_pos[N];
	double x;
	double y;


	// Set all init pos to 1
	memset(init_x_pos,-1.0,N*sizeof(double));
	memset(init_y_pos,-1.0,N*sizeof(double));
	for(int i = 0; i < N ; i++){
		// unique id
		p_array[i].id = i;

		// position init
		x = rand_pos_component();
		y = rand_pos_component();
		while(checkArray(init_x_pos,x,i) == 1 ){
			x = rand_pos_component();
		}
		while(checkArray(init_y_pos,y,i) == 1 ){
			y = rand_pos_component();
		}
		p_array[i].pos_x = x;
		p_array[i].pos_y = y;
		// Save positions in our list to be checked later
		init_x_pos[i] = x;
		init_y_pos[i] = y;

		// vel and acc init
		p_array[i].vel_x = rand_val_component();
		p_array[i].vel_y = rand_val_component();
		p_array[i].acc_x = rand_val_component();
		p_array[i].acc_y = rand_val_component();
		p_array[i].acc_x_next_step = p_array[i].acc_x;
		p_array[i].acc_y_next_step = p_array[i].acc_y;

	}
}


/*
	The following method intializes particals uniformly inside L_n X L_n box 
	** UPTO 100 particles max using this method **
	-------------
	| *	* *	* *	|
	| *	* * * *	|
	| * * * * * |
	-------------
	*/
void initialize_particles_uniform(struct Particle p_array[N]){
	
	double padding = 1.0;
	double p_dist = (L_d) / ceil(sqrt(N)); // Divide up particles


	// Generate all possible coordinate and store in this list
	int num_possible_coordinates = (int)pow(ceil(sqrt(N)),2);
	double x_list[num_possible_coordinates];
	double y_list[num_possible_coordinates];

	int idx = 0;
		for(double x = padding/2.0; x <= (L_d -(padding/2.0)); x += p_dist ){
			for(double y = padding/2.0; y <= (L_d -(padding/2.0)); y += p_dist ){

					p_array[idx].id = idx;
					p_array[idx].pos_x = x;
					p_array[idx].pos_y = y;

					// Init vel acc
					p_array[idx].vel_x = rand_val_component();
					p_array[idx].vel_y = rand_val_component();
					p_array[idx].acc_x = rand_val_component();
					p_array[idx].acc_y = rand_val_component();
					p_array[idx].acc_x_next_step = p_array[idx].acc_x;
					p_array[idx].acc_y_next_step = p_array[idx].acc_y;
					idx++;

					if(idx >= N){

						return;
					}
			}
		}
}

// Entry Point into program
int main(){
	srand(time(NULL)); // Random number whose seed is dependent on time
	clock_t begin = clock(); // count cpu time

	struct Particle p_array[N]; // Hold all particle structs
	initialize_particles_uniform(p_array); // Choose a particle intialization method

	char Filename[21] = FILENAME; //Filename

	/*
	Use following lines to generate table of particles and associated init data

	printf("\nParticle Table at t = 0\n\n");
	particlePrinter(p_array);
	printf("\n");
	*/

	// Compute statisitcs
	double init_int_energy = Internal_energy(p_array);
	double rolling_avg_vel = 0.0;
	double rolling_avg_sep_dist = 0.0;
	double velocity_stat[SIM_TIME]; // Store v^2

	clock_t end = clock(); // count cpu time
	printf("\nInitialization Time     : %e seconds\n", (double)(end-begin) / CLOCKS_PER_SEC );

	int miss_cnt;


	for(int t = 0; t < SIM_TIME; t++){
		UpdatePosition(p_array);
		UpdateAcceleration(p_array);
		UpdateVelocity(p_array);
		enforce_boundary(p_array);
		// Print current time in place
		printf("\rCurrent Time: %d timesteps",t);
		fflush(stdout);
		miss_cnt = misbehaving_particle_counter(p_array);
		


		// collect stats
		rolling_avg_vel+= Avg_velocity(p_array);
		rolling_avg_sep_dist += Avg_Seperation(p_array);

		velocity_stat[t] = pow(Avg_velocity(p_array),2.0);
	}

	//averaging stats
	rolling_avg_vel /= SIM_TIME;
	rolling_avg_sep_dist /= SIM_TIME;

	double final_int_energy = Internal_energy(p_array);

	WriteToFile(Filename,SIM_TIME,velocity_stat); // Write file containing v^2 vs time

	end = clock(); // count cpu time


	printf("\r"); // remove current time
	printf("Program Execution Time  : %e seconds\n", (double)(end-begin) / CLOCKS_PER_SEC );
	printf("TOTAL SIMULATION TIME   : %e\n",dt*SIM_TIME );
	printf("Avg_Seperation Distance : %f\n",rolling_avg_sep_dist );
	printf("Internal Energy         : %f\n", final_int_energy );
	printf("Average Particle Speed  : %f\n",rolling_avg_vel );

	

	// Check for errors and create a report
	printf("\nERROR REPORT:\n");
	printf("PARTICLES WITHIN BOX : ");
	if(miss_cnt > 0){
			printf("[FAIL] %d Particles outside box\n",miss_cnt);
	} else {
		printf("[OK]\n");
	}
	printf("ENERGY CONSERVATION  : ");
	if( (final_int_energy/init_int_energy) <= 1.5 && (final_int_energy/init_int_energy) >= 0.5){
		printf("[OK]\n" );
	} else {
		printf("[FAIL] Initial energy: %f Final energy %f \n",init_int_energy,final_int_energy );
	}

	printf("\nData Output Saved As %s\n",Filename);


	return 1;
}
