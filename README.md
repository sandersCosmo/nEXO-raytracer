# nEXO-raytracer
A ray tracer for the SLAC nEXO group's Liquid Xenon Cell

## Running
The script requires you to define the number of decays you want to simulate as an additional argument when executing the script, as well as an ID to be used in naming the data file. The format for this execution is ```python nEXO_raytracer <n_decays> <file_id>```. For example, if ```n_decays=100``` and ```file_id=0```, then the counts for 100 decays in the bottom detector are written to ```det1_counts0.npy```.

The data files are written to the ```data/``` directory in the same directory as the script. Either create this directory or change the location of the saved files in the ```np.save``` statements at the end of the ```loop``` function.

### Using SDF
SDF is the shared computing resource at SLAC. It uses the job manager slurm. There is an example bash file that is used to schedule a job. To submit the job, execute ```sbatch <jobfile>```. For running many jobs, I'd recommend making a script that creates many of these job files and executes the ```sbatch``` command. The detailed slurm documentation is located here: https://slurm.schedmd.com/.

## Geometry
Some operations need derivations, which are given here.

### Reflection off Cylinder
All reflections in the code are diffuse, so when picking a new random angle for the ray to be reflected in, we need to be sure we only choose angles that stay inside the cylinder. This is to say, only the half of angles on the inside of the line tangent to the point on the circle.

![IMG_3188](https://user-images.githubusercontent.com/110847459/183491339-abb9aeb2-7ed3-4866-8899-50d7402e4de9.jpg)

It is important to note that ```arctan``` returns values along $[-\pi/2, \pi/2]$. If $x<0$, then we must add $\pi$ to $\phi$ in order to get the correct $\phi$.

### Refraction
The rays refract through the liquid/gas boundary. So, we just work out the law of refraction in an arbitrary direction in 3 dimensions. (note that "Snell's Law" is misnamed, it was discovered 100s of years ealrier by Persian scientist Ibn Sahl.)

![IMG_3189](https://user-images.githubusercontent.com/110847459/183511921-3af464dd-2726-432b-a499-4fe4ab10ff41.jpg)

Here we always demand that the magnitude of the velocity vector be 1. Similar to before, ```arctan``` returns values along $[-\pi/2, \pi/2]$. If the incorrect angle is chosen, it will only result in a sign error. So, the code will flip the sign to the correct value (ie. refraction will never change the sign of a component.)

It is possible that the argument to the ```arcsin``` in $dz$ is greater than 1. If this is the case, total internal reflection takes place, so the sign of $dz$ is simply flipped.

## References
PTFE Reflectivity in Liquid Xenon: https://doi.org/10.1063/1.3318681
Liquid Xenon Attenuation: https://arxiv.org/abs/physics/0307044#:~:text=The%20attenuation%20length%20and%20refractive,to%20be%201.69%20%2B%2D%200.02.
