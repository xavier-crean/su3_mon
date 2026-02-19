import numpy as np
import itertools
from gtda.homology import CubicalPersistence
from opt_einsum import contract
from multiprocessing import Pool
from sys import argv
from configuration_u1 import *
from multiprocessing import Pool
import os

"""
On command line: 
argv[1] = Nt.Nx.Ny.Nz, argv[2] = beta, argv[3] = N_meas, argv[4] = sweep_step, argv[5] = sweep_start, argv[6] = N_processors,
E.g., "python3 ./src/betti_compute.py 8.64.64.64 5.8900 400 2000 10000 5"
"""

PLANE_LABELS = {(0, 1): 0, (0, 2): 1, (0, 3): 2, (1, 2): 3, (1, 3): 4, (2, 3): 5}
NCOLOR = 3

def diracStringNum(plaq_ang):
    """Take a plaquette angle and return the number of Dirac strings passing through as int"""
    if (-4 * np.pi) < plaq_ang <= (-3 * np.pi):
        return -2
    elif (-3 * np.pi) < plaq_ang <= (-1 * np.pi):
        return -1
    elif (-1 * np.pi) < plaq_ang <= np.pi:
        return 0
    elif np.pi < plaq_ang <= (3 * np.pi):
        return 1
    elif (3 * np.pi) < plaq_ang <= (4 * np.pi):
        return 2
    else:
        raise Exception("Plaq value too large/small.")


def levi_cevita(dim):
    """Generates Levi-Cevita symbol for dimensions dim"""
    arr = np.zeros(tuple([dim for _ in range(dim)]))
    for x in itertools.permutations(tuple(range(dim))):
        mat = np.zeros((dim, dim), dtype=np.int32)
        for i, j in zip(range(dim), x):
            mat[i, j] = 1
        arr[x] = int(np.linalg.det(mat))
    return arr


# Intialise 4d levi-cevita tensor
levi = levi_cevita(4)


def mon(plang):
    """
    Input: array (dtype=float64) of plaquette angles
    Return: array (dtype=int32) of monopole currents on dual lattice
    """
    shape = plang.shape[:4]
    
    # Dirac string array
    n = np.zeros(shape + (4, 4), dtype=np.int32)
    for t, x, y, z in np.ndindex(shape):
        for k, v in PLANE_LABELS.items():
            n[t, x, y, z, k[0], k[1]] = diracStringNum(plang[t, x, y, z, v])

    # monopole 3-cube array
    M = np.zeros(shape + (4,), dtype=np.int32)
    for t, x, y, z in np.ndindex(shape):
        # free index
        for r in range(4):
            total = 0
            for s in range(4):
                if s == 0:
                    ds = [1, 0, 0, 0]
                elif s == 1:
                    ds = [0, 1, 0, 0]
                elif s == 2:
                    ds = [0, 0, 1, 0]
                elif s == 3:
                    ds = [0, 0, 0, 1]
                n_diff = (
                    n[
                        (t + ds[0]) % shape[0],
                        (x + ds[1]) % shape[1],
                        (y + ds[2]) % shape[2],
                        (z + ds[3]) % shape[3],
                    ]
                    - n[t, x, y, z]
                )
                total += contract("mn,mn->", levi[r, s], n_diff)
            M[t, x, y, z][r] = total

    # monopole dual 1-cube array
    j = np.zeros(shape + (4,), dtype=np.int32)
    for t, x, y, z in np.ndindex(shape):
        for r in range(4):
            if r == 0:
                dr = [1, 0, 0, 0]
            elif r == 1:
                dr = [0, 1, 0, 0]
            elif r == 2:
                dr = [0, 0, 1, 0]
            elif r == 3:
                dr = [0, 0, 0, 1]
            j[t, x, y, z][r] = M[
                (t + dr[0]) % shape[0],
                (x + dr[1]) % shape[1],
                (y + dr[2]) % shape[2],
                (z + dr[3]) % shape[3],
            ][r]
    return j


# +- unit vector in each direction
pm = np.array(
    [
        [[1, 0, 0, 0], [-1, 0, 0, 0]],
        [[0, 1, 0, 0], [0, -1, 0, 0]],
        [[0, 0, 1, 0], [0, 0, -1, 0]],
        [[0, 0, 0, 1], [0, 0, 0, -1]],
    ]
)


def cubicalFiltration(m):
    """
    Trivial filtration:
        s=0: non-zero current lines and their boundary vertices
        s=np.inf: all other d-cubes enter completing the 4-torus
    """
    filt = np.zeros(tuple([2 * s for s in m.shape[:4]]))
    shape = filt.shape

    for t, x, y, z in np.ndindex(shape):
        dim = t % 2 + x % 2 + y % 2 + z % 2
        if dim == 1:
            parDirs = tuple([i for i in range(4) if [t, x, y, z][i] % 2 == 1])
            j = m[
                (t - (t % 2)) // 2,
                (x - (x % 2)) // 2,
                (y - (y % 2)) // 2,
                (z - (z % 2)) // 2,
                parDirs[0],
            ]
            if j != 0:
                filt[t, x, y, z] = 0.0
            else:
                filt[t, x, y, z] = np.inf
        elif dim == 2:
            filt[t,x,y,z] = np.inf
        elif dim == 3:
            filt[t,x,y,z] = np.inf
        elif dim == 4:
            filt[t,x,y,z] = np.inf

    for t,x,y,z in np.ndindex(shape):
        if (t%2 + x%2 + y%2 + z%2) == 0:
            dirs = np.concatenate(pm[[i for i in range(4) if [t,x,y,z][i]%2 == 0]]) # 0 since looking up dimension
            filt[t,x,y,z] = min([filt[(t+dt)%shape[0],(x+dx)%shape[1],(y+dy)%shape[2],(z+dz)%shape[3]] for [dt,dx,dy,dz] in dirs])

    return filt


def process_comb(confs_0, confs_1):
    return np.abs(mon(plaquettes(confs_0))) + np.abs(mon(plaquettes(confs_1)))


def multiplicity(diagram, homology_dimensions=None):
    if homology_dimensions is None:
        homology_dimensions = np.unique(diagram[:, 2])  # Extract unique dimensions
    counts = []
    for dim in homology_dimensions:
        subdiagram = diagram[diagram[:, 2] == dim]
        subdiagram = subdiagram[subdiagram[:,0] != subdiagram[:,1]]
        counts.append(len(subdiagram))
    return np.array(counts)

### Initialise cubical persistence object
cp = CubicalPersistence(
            homology_dimensions=[0, 1],
            coeff=2,
            periodic_dimensions=np.array([True, True, True, True]),
            n_jobs=1  # Local computation per worker
        )

def process_config_pair(task):
    """Pipeline for a single configuration pair to compute Betti numbers -- with explicit parameters"""
    # Unpack parameters
    config_idx, beta, Nt, Nx, Ny, Nz, start, sweep_step = task
    
    try:
        # Calculate exact file index using passed parameters
        r = start + config_idx * sweep_step
        
        # Read files using path parameters
        fn0 = f"data/cnfgs/{Nt}.{Nx}.{Ny}.{Nz}/{beta:.4f}/conf_u1_subg0.dat{r}"
        fn1 = f"data/cnfgs/{Nt}.{Nx}.{Ny}.{Nz}/{beta:.4f}/conf_u1_subg1.dat{r}"
        
        with open(fn0, "rb") as f:
            conf0 = parseConfig(f.read(), Nt, Nx, Ny, Nz)
        with open(fn1, "rb") as f:
            conf1 = parseConfig(f.read(), Nt, Nx, Ny, Nz)
        
        # 2. Process data
        mon_comb = process_comb(conf0, conf1)
        filt = cubicalFiltration(mon_comb)
        
        # 3. Compute persistence
        pd = cp.fit_transform([filt])[0]
        
        # 4. Compute Betti numbers
        betti = multiplicity(pd)
        
        # Return index with results
        return (config_idx, betti)
        
    except Exception as e:
        print(f"Error processing config {config_idx}: {str(e)}")
        return (config_idx, np.zeros(2, dtype=int))


def main():
    ### Parameter set-up
    param = np.array(argv[1].split("."), dtype=int)
    Nt, Nx, Ny, Nz = param[0], param[1], param[2], param[3]
    b = float(argv[2])
    N = int(argv[3])
    sweep_step = int(argv[4])
    start = int(argv[5])
    NP = int(argv[6])
    
    ### File paths
    os.makedirs(f"../data/observables/betti/{Nt}.{Nx}.{Ny}.{Nz}/", exist_ok=True)
    final_filename = f"../data/observables/betti/{Nt}.{Nx}.{Ny}.{Nz}/betti_NtNsNsNs={Nt}{Nx}{Ny}{Nz}_b={b:.4f}.csv"
    temp_filename = f"../data/observables/betti/{Nt}.{Nx}.{Ny}.{Nz}/betti_temp_NtNsNsNs={Nt}{Nx}{Ny}{Nz}_b={b:.4f}.csv"
    
    # Check if final output already exists
    if os.path.exists(final_filename):
        print(f"Final output {final_filename} already exists. Exiting.")
        return
    
    ### Initialize or load progress
    if os.path.exists(temp_filename):
        try:
            h = np.loadtxt(temp_filename, delimiter=",", dtype=int)
            if h.shape != (N, 2):
                raise ValueError("Temp file has incorrect shape")
            print(f"Resuming from existing temp file {temp_filename}")
        except:
            print("Temp file corrupt/invalid. Starting fresh.")
            h = np.zeros((N, 2), dtype=int)
    else:
        h = np.zeros((N, 2), dtype=int)
    
    ### Find remaining configurations to process
    remaining_indices = [i for i in range(N) if np.all(h[i] == 0)]
    if not remaining_indices:
        print("All configurations already processed. Moving temp to final.")
        np.savetxt(final_filename, h, delimiter=",", fmt="%d")
        if os.path.exists(temp_filename):
            os.remove(temp_filename)
        return
    
    ### Create task list for remaining configurations
    tasks = [
        (i, b, Nt, Nx, Ny, Nz, start, sweep_step)
        for i in remaining_indices
    ]
    
    ### Process configurations with periodic saving
    processed_count = 0
    total_to_process = len(remaining_indices)
    
    with Pool(processes=NP) as pool:
        # Process tasks and handle results with index mapping
        results = pool.imap(process_config_pair, tasks, chunksize=1)
        
        try:
            for i_global, betti in results:
                h[i_global] = betti
                processed_count += 1
                
                # Save progress every NP processed configurations
                if processed_count % NP == 0:
                    np.savetxt(temp_filename, h, delimiter=",", fmt="%d")
                    print(f"Checkpoint: {processed_count}/{total_to_process} processed (Total completed: {N - len(remaining_indices) + processed_count}/{N})")
        
        except Exception as e:
            print(f"Error during processing: {str(e)}")
            print("Saving current progress before exiting...")
            np.savetxt(temp_filename, h, delimiter=",", fmt="%d")
            raise
    
    ### Final save and cleanup
    np.savetxt(final_filename, h, delimiter=",", fmt="%d")
    if os.path.exists(temp_filename):
        os.remove(temp_filename)
    print(f"Completed all {N} configurations. Results saved to {final_filename}")

if __name__ == "__main__":
    main()