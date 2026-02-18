from struct import unpack
import itertools
import numpy as np

PLANE_LABELS = {(0, 1): 0, (0, 2): 1, (0, 3): 2, (1, 2): 3, (1, 3): 4, (2, 3): 5}


def parseConfig(config_buffer, Nt, Nx, Ny, Nz):
    DTYPE_LEN = 8
    STDIM, DSIZE_T, DSIZE_X, DSIZE_Y, DSIZE_Z = 4, Nt, Nx, Ny, Nz

    def AngleFromBytes(buffer):
        return unpack("<d", buffer)[0]

    def cart_to_lex(t, x, y, z):
        return (
            t
            + (x * DSIZE_T)
            + (y * DSIZE_T * DSIZE_X)
            + (z * DSIZE_T * DSIZE_X * DSIZE_Y)
        )

    config_array = np.array(
        [
            [
                [
                    [
                        [
                            AngleFromBytes(
                                config_buffer[
                                    DTYPE_LEN
                                    * (
                                        STDIM * cart_to_lex(t, x, y, z) + link
                                    ) : DTYPE_LEN
                                    * (STDIM * cart_to_lex(t, x, y, z) + link + 1)
                                ]
                            )
                            for link in range(STDIM)
                        ]
                        for z in range(DSIZE_Z)
                    ]
                    for y in range(DSIZE_Y)
                ]
                for x in range(DSIZE_X)
            ]
            for t in range(DSIZE_T)
        ]
    )
    return config_array


def wilsonLoop(confang, t, x, y, z, directions):
    d_mu = [[1 if i == d else 0 for i in range(4)] for d in directions]

    ang = (
        confang[t, x, y, z, directions[0]]
        + confang[
            (t + d_mu[0][0]) % confang.shape[0],
            (x + d_mu[0][1]) % confang.shape[1],
            (y + d_mu[0][2]) % confang.shape[2],
            (z + d_mu[0][3]) % confang.shape[3],
            directions[1],
        ]
        - confang[
            (t + d_mu[1][0]) % confang.shape[0],
            (x + d_mu[1][1]) % confang.shape[1],
            (y + d_mu[1][2]) % confang.shape[2],
            (z + d_mu[1][3]) % confang.shape[3],
            directions[0],
        ]
        - confang[t, x, y, z, directions[1]]
    )
    return np.array(ang)


def plaquettes(conf):
    plaq = np.zeros(conf.shape[:4] + (6,), dtype=np.double)
    for t, x, y, z in itertools.product(*[range(s) for s in conf.shape[:4]]):
        for k, v in PLANE_LABELS.items():
            plaq[t, x, y, z, v] = wilsonLoop(conf, t, x, y, z, k)
    return plaq
