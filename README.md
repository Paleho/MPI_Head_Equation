# Thermal Equation - MPI 

### Contributors
#### - Poutas Sokratis (poutasok@gmail.com)
#### - Tzomaka Afroditi (afrodititzomaka@gmail.com)

---

### Here are three parallelized algorithms that solve the heat equation using **OpenMPI**

- Jacobi
- Gauss Seidel 
- Red-Black

### All three implementations are designed for *distributed memory parallel systems*

In every directory there is our source code and a Makefile that produces 2 executables: one that terminates after a convergence check is met and one that runs for 256 iterations

> run `make` in every directory to compile the code

to run the code:

```
module load openmpi/1.8.3
mpirun -np 64 --mca btl self,tcp jacobi_conv 1024 1024 8 8
``` 
(for running jacobi with 64 processes in a 1024x1024 table with 8x8 process grid)


## Preformance Notes

- Jacobi is the most parallelized version and this makes it perform better in fixed iterations scenarios
- Gauss Seidel utilizes less parallelization but converges way faster
- Red-Black converges even faster but has the worse performance in fixed iterations scenarios