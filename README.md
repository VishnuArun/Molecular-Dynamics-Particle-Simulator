# Molecular Dynamics Particle Simulator

This simulation will compute the internal energy, average interparticle seperation, and average velocity of particles in a 2D box of specified size.



- Particles interact via a standard Lennard-Jones potential. 
- $$ V(r) = 4 \epsilon \left( \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r}\right)^{6} \right) $$


- Velocity-Verlet Algorithm used for time integration
- $$ x(t+\Delta t) = xt) + v(t) + \frac{1}{2}a(t)\Delta t^2 $$
- $$ v(t + \Delta t) = v(t) + \frac{a(t) + a(t + \Delta t)}{2}\Delta t$$

- System held a constant temperature T. Accomplished by applying boundary conditions When a particle collides with any of the walls, the tangential component of the velocity of the particle is conserved, and the incident normal velocity is updated into an outgoing normal velocity given by [ [van Beijeren](http://arxiv.org/abs/1411.2983), 2014]. 
- $$v_{n}^{\prime} = \sqrt{ - \frac{2}{\beta m} \ln \left( 1 - {\rm exp} \left( - \frac{\beta m v_{n}^{2}}{2} \right) \right) }$$ 

Mass of the particles are all assumed identical. Dimensions are transformed into unitless quantities.


### Considerations regarding implementation.
- Intialization of particles: We have two seperate modules for particle intialization. 
	- Method 1: The particles are randomly assigned and ensure they do not overlap. This can cause particles to be very close to one another especially when simulating a greater number of particles (N > 30) in a relatively small box (L < 15).
	- Method 2: The particles are uniformly spaced. This implementation allows the simulation to run on a larger number of particles, but is unrealistic as a physical condition.

- Boundary conditions with larger dt time increments. The dt time increment is the smallest amount of time that can pass between timesteps. If this variable is large the particles can develop large velcoties and leave the box very quickly. We have implemented functions that will correct particles to stay contained within the box. They "reflect" to a distance inside the box equivalent to the distance of the particle from a wall.

### Error Reporting Mechanisms
We have implemented a handful of Error Reporting Mechanisms to ensure the simulation follows real-world physical limitations. All errors will be reported after completion of the simulation.

- Conservation of Energy: One of the key assumption we have held is that the temperature is constant, that is, no energy can be lost resulting in an impact of a particle on wall (Perfectly elastic collision). We have developed a module to monitor the total system energy and report whether it has changed significantly (By an order of magniture) from the intial internal energy. The internal energy is computed as the kinetic + potential energy of all particles.
- Particles cannot "walk through walls" that is have position coordinates outside the dimensions of the box. We ensure every position update for all particles requires this condition be met before updating. If the particle is found to exit the box as a result of a future update we modify the position to act as a "refletion" off the nearest boundary wall. Boundary conditions are therby enforced at every position update. A module will monitor this behavior and report if a particle has gone out of bounds.


### Statistic Report
Once the simulation time has expired the program will produce a report. This report contains the simulation time, average seperation distance, Internal energy, and average particle speed. The program will also output a .dat file containing the average velocity squared of all particles eveolved over time. 


