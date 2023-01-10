import numpy as np
from numpy import zeros
import math


def makeGrid(
    nRows,
    nTurbines,
    dx,
    dy,
    shear,
    rotate,
    turbs_per_row,
    x_start,
    y0,
    turbineX,
    turbineY,
):

    nRows = int(nRows)
    nTurbines = int(nTurbines)
    dx = float(dx)
    dy = float(dy)
    shear = float(shear)
    rotate = float(rotate)
    turbs_per_row = np.array(turbs_per_row, dtype=int)
    x_start = np.array(x_start, dtype=float)
    y0 = float(y0)
    turbineX = np.zeros(nTurbines, dtype=float)
    turbineY = np.zeros(nTurbines, dtype=float)

    # local
    index = 0
    i_fl = 1.0
    x = np.zeros(nTurbines, dtype=float)
    y = np.zeros(nTurbines, dtype=float)
    xc = np.zeros(nTurbines, dtype=float)
    yc = np.zeros(nTurbines, dtype=float)

    for i in range(nRows):
        j_fl = 1.0
        for j in range(turbs_per_row[i]):
            x[index] = x_start[i] * dx + dx * j_fl + i_fl * shear
            y[index] = y0 + dy * i_fl
            index += 1
            j_fl += 1.0
        i_fl += 1.0

    minx = 100000000.0
    maxx = -100000000.0
    miny = 100000000.0
    maxy = -100000000.0
    for i in range(nTurbines):
        if x[i] < minx:
            minx = x[i]
        if x[i] > maxx:
            maxx = x[i]
        if y[i] < miny:
            miny = y[i]
        if y[i] > maxy:
            maxy = y[i]

    xc = x - (maxx + minx) / 2
    yc = y - (maxy + miny) / 2
    rotate_rad = (rotate * np.pi) / 180.0
    turbineX = np.cos(rotate_rad) * xc - np.sin(rotate_rad) * yc
    turbineY = np.sin(rotate_rad) * xc + np.cos(rotate_rad) * yc

    return turbineX, turbineY


def makegrid_fortran_dv(
    nrows,
    nturbines,
    dx,
    dxd,
    dy,
    dyd,
    shear,
    sheard,
    rotate,
    rotated,
    turbs_per_row,
    x_start,
    y0,
    turbinex,
    turbinexd,
    turbiney,
    turbineyd,
    nbdirs,
):
    #  Hint: nbdirs should be the maximum number of differentiation directions
    # in
    # out
    # local
    index = 1
    i_fl = 1.0
    xd = np.zeros((nbdirs, nturbines), dtype=np.float64)
    yd = np.zeros((nbdirs, nturbines), dtype=np.float64)
    xcd = np.zeros((nbdirs, nturbines), dtype=np.float64)
    ycd = np.zeros((nbdirs, nturbines), dtype=np.float64)
    x = np.zeros(nturbines, dtype=float)
    y = np.zeros(nturbines, dtype=float)
    xc = np.zeros(nturbines, dtype=float)
    yc = np.zeros(nturbines, dtype=float)

    for i in range(1, nrows + 1):
        j_fl = 1.0
        for j in range(1, turbs_per_row[i - 1] + 1):
            for nd in range(1, nbdirs + 1):
                xd[nd - 1, index - 1] = (
                    x_start[i - 1] * dxd[nd - 1]
                    + j_fl * dxd[nd - 1]
                    + i_fl * sheard[nd - 1]
                )
                yd[nd - 1, index - 1] = i_fl * dyd[nd - 1]
            x[index - 1] = x_start[i - 1] * dx + dx * j_fl + i_fl * shear
            y[index - 1] = y0 + dy * i_fl
            index = index + 1

    minx = 100000000.0
    maxx = -100000000.0
    miny = 100000000.0
    maxy = -100000000.0
    minxd = zeros(nbdirs)
    maxxd = zeros(nbdirs)
    minyd = zeros(nbdirs)
    maxyd = zeros(nbdirs)
    rotate_radd = zeros(nbdirs)
    for i in range(1, nturbines):
        if x[i] < minx:
            for nd in range(1, nbdirs):
                minxd[nd] = xd[nd, i]
            minx = x[i]
        if x[i] > maxx:
            for nd in range(1, nbdirs):
                maxxd[nd] = xd[nd, i]
            maxx = x[i]
        if y[i] < miny:
            for nd in range(1, nbdirs):
                minyd[nd] = yd[nd, i]
            miny = y[i]
        if y[i] > maxy:
            for nd in range(1, nbdirs):
                maxyd[nd] = yd[nd, i]
            maxy = y[i]

    xc = x - (maxx + minx) / 2
    yc = y - (maxy + miny) / 2
    rotate_rad = rotate * math.pi / 180.0
    for nd in range(1, nbdirs):
        xcd[nd, :] = xd[nd, :] - (maxxd[nd] + minxd[nd]) / 2
        ycd[nd, :] = yd[nd, :] - (maxyd[nd] + minyd[nd]) / 2
        rotate_radd[nd] = math.pi * rotated[nd] / 180.0
        turbinexd[nd, :] = (
            np.cos(rotate_rad) * xcd[nd, :]
            - rotate_radd[nd] * np.sin(rotate_rad) * xc
            - rotate_radd[nd] * np.cos(rotate_rad) * yc
            - np.sin(rotate_rad) * ycd[nd, :]
        )
        turbineyd[nd, :] = (
            rotate_radd[nd] * np.cos(rotate_rad) * xc
            + np.sin(rotate_rad) * xcd[nd, :]
            + np.cos(rotate_rad) * ycd[nd, :]
            - rotate_radd(nd) * np.sin(rotate_rad) * yc
        )

    turbinex = np.cos(rotate_rad) * xc - np.sin(rotate_rad) * yc
    turbiney = np.sin(rotate_rad) * xc + np.cos(rotate_rad) * yc

    return turbinex, turbiney, turbinexd, turbineyd
