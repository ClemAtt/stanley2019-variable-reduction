import math
import numpy as np
from numpy import zeros


def makegrid_fortran(nRows, nGrid, dx, dy, shear, rotate, turbs_per_row, x_start, y0):
    """Makes a grid of turbines

    Args:
        nRows (int): Number of rows of turbines
        nGrid (int): Number of turbines in the grid
        dx (float): Change in x direction
        dy (float): Change in y direction
        shear (float): Shear of grid
        rotate (float): Rotation of grid
        turbs_per_row (int): Turbines per row
        x_start (float): Starting x position
        y0 (float): Starting y position

    Returns:
        Turbine x and y positions: Array of turbine x and y positions
    """
    # print("nRows: ", nRows)
    # print("nGrid: ", nGrid)
    # print("turbs_per_row: ", turbs_per_row)
    turbineX, turbineY = [], []
    i_fl = 1.0
    x = np.zeros(nGrid, dtype=np.float64)
    y = np.zeros(nGrid, dtype=np.float64)

    for i in range(nRows):
        j_fl = 1.0
        for j in range(turbs_per_row[i]):
            x[i + j] = x_start[i] * dx + dx * j_fl + i_fl * shear
            y[i + j] = y0 + dy * i_fl
            j_fl += 1.0
        i_fl += 1.0

    minx, maxx, miny, maxy = 1e8, -1e8, 1e8, -1e8
    for i in range(nGrid):
        minx = min(minx, x[i])
        maxx = max(maxx, x[i])
        miny = min(miny, y[i])
        maxy = max(maxy, y[i])

    xc = x - (maxx + minx) / 2
    yc = y - (maxy + miny) / 2
    rotate_rad = (rotate * math.pi) / 180.0
    turbineX = math.cos(rotate_rad) * xc - math.sin(rotate_rad) * yc
    turbineY = math.sin(rotate_rad) * xc + math.sin(rotate_rad) * yc
    return turbineX, turbineY


def makeBoundary(nBounds, nOuter, boundX, boundY, start):
    turbineX, turbineY = [], []
    lenBound = [0] * (nBounds - 1)
    for i in range(nBounds - 1):
        lenBound[i] = math.sqrt(
            (boundX[i + 1] - boundX[i]) ** 2 + (boundY[i + 1] - boundY[i]) ** 2
        )

    circumference = sum(lenBound)
    spacing = circumference / nOuter
    bound_loc = start

    for i in range(nOuter):
        while bound_loc > circumference:
            bound_loc -= circumference
        while bound_loc < 0.0:
            bound_loc += circumference
        done = 0
        for j in range(nBounds):
            if done == 0:
                if bound_loc < sum(lenBound[:j]):
                    if j == 1:
                        x1, y1, x2, y2 = (
                            boundX[nBounds - 1],
                            boundY[nBounds - 1],
                            boundX[j],
                            boundY[j],
                        )
                    else:
                        x1, y1, x2, y2 = (
                            boundX[j - 2],
                            boundY[j - 2],
                            boundX[j - 1],
                            boundY[j - 1],
                        )
                    slope = (y2 - y1) / (x2 - x1)
                    intersect = y1 - x1 * slope
                    turbineX.append(bound_loc / slope + intersect)
                    turbineY.append(bound_loc * slope + intersect)
                    done = 1
                    bound_loc += spacing
    return turbineX, turbineY


def turbineLocations(
    nBounds,
    nRows,
    nTurbines,
    nOuter,
    nGrid,
    dx,
    dy,
    shear,
    rotate,
    turbs_per_row,
    x_start,
    y0,
    start,
    boundX,
    boundY,
):
    outerX, outerY = makeBoundary(nBounds, nOuter, boundX, boundY, start)
    innerX, innerY = makegrid_fortran(
        nRows, nGrid, dx, dy, shear, rotate, turbs_per_row, x_start, y0
    )
    turbineX, turbineY = [], []
    turbineX[:nOuter] = outerX
    turbineY[:nOuter] = outerY
    turbineX[nOuter:] = innerX
    turbineY[nOuter:] = innerY
    return turbineX, turbineY


def makegrid_fortran_dv(
    nrows,
    ngrid,
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
    x = np.zeros(ngrid, dtype=np.float64)
    y = np.zeros(ngrid, dtype=np.float64)

    index = 1
    i_fl = 1.0
    xd = np.zeros((nbdirs, ngrid), dtype=np.float64)
    yd = np.zeros((nbdirs, ngrid), dtype=np.float64)

    for i in range(nrows):
        j_fl = 1.0
        for j in range(turbs_per_row[i]):
            for nd in range(nbdirs):
                xd[nd][index] = (
                    x_start[i] * dxd[nd] + j_fl * dxd[nd] + i_fl * sheard[nd]
                )
                yd[nd][index] = i_fl * dyd[nd]
            x[index] = x_start[i] * dx + dx * j_fl + i_fl * shear
            y[index] = y0 + dy * i_fl
            index = index + 1
            j_fl = j_fl + 1.0
        i_fl = i_fl + 1.0

    minx = 100000000.0
    maxx = -100000000.0
    miny = 100000000.0
    maxy = -100000000.0
    minxd = np.zeros(nbdirs, dtype=np.float64)
    minyd = np.zeros(nbdirs, dtype=np.float64)
    maxxd = np.zeros(nbdirs, dtype=np.float64)
    maxyd = np.zeros(nbdirs, dtype=np.float64)
    for i in range(ngrid):
        if x[i] < minx:
            minxd = np.copy(xd[:, i])
            minx = x[i]
        if x[i] > maxx:
            maxxd = np.copy(xd[:, i])
            maxx = x[i]
        if y[i] < miny:
            minyd = np.copy(yd[:, i])
            miny = y[i]
        if y[i] > maxy:
            maxyd = np.copy(yd[:, i])
            maxy = y[i]

    rotate_rad = rotate * np.pi / 180.0
    rotate_radd = rotated * np.pi / 180.0

    xc = np.copy(x)
    yc = np.copy(y)

    xc = x - (maxx + minx) / 2
    yc = y - (maxy + miny) / 2
    rotate_rad = rotate * 3.1415926535 / 180.0
    for nd in range(nbdirs):
        xcd = xd[nd, :] - (maxxd[nd] + minxd[nd]) / 2
        ycd = yd[nd, :] - (maxyd[nd] + minyd[nd]) / 2
        rotate_radd = 3.1415926535 * rotated[nd] / 180.0
        turbinexd[nd, :] = (
            np.cos(rotate_rad) * xcd
            - rotate_radd * np.sin(rotate_rad) * xc
            - rotate_radd * np.cos(rotate_rad) * yc
            - np.sin(rotate_rad) * ycd
        )
        turbineyd[nd, :] = (
            rotate_radd * np.cos(rotate_rad) * xc
            + np.sin(rotate_rad) * xcd
            + np.cos(rotate_rad) * ycd
            - rotate_radd * np.sin(rotate_rad) * yc
        )
    turbinex = np.cos(rotate_rad) * xc - np.sin(rotate_rad) * yc
    turbiney = np.sin(rotate_rad) * xc + np.cos(rotate_rad) * yc

    return turbinex, turbinexd, turbiney, turbineyd


def makeboundary_dv(
    nbounds,
    nouter,
    boundx,
    boundy,
    start,
    startd,
    turbinex,
    turbinexd,
    turbiney,
    turbineyd,
    nbdirs,
):
    #  Hint: nbdirs should be the maximum number of differentiation directions
    for i in range(1, nbounds - 1 + 1):
        arg1 = (boundx(i + 1) - boundx(i)) ** 2 + (boundy(i + 1) - boundy(i)) ** 2
        lenbound[i] = sqrt(arg1)
    circumference = sum(lenbound)
    spacing = circumference / nouter
    for nd in range(1, nbdirs + 1):
        bound_locd[nd] = startd[nd]
    bound_loc = start
    turbinexd[:, :] = 0.0_8
    turbineyd[:, :] = 0.0_8
    for i in range(1, nouter + 1):
        while bound_loc > circumference:
            bound_loc = bound_loc - circumference
        while bound_loc < 0.0:
            bound_loc = bound_loc + circumference
        done = 0
        for j in range(1, nbounds + 1):
            if done == 0:
                if bound_loc < sum(lenbound[:j]):
                    for nd in range(1, nbdirs + 1):
                        turbinexd[nd, i] = (
                            (boundx(j + 1) - boundx(j)) * bound_locd[nd] / lenbound[j]
                        )
                        turbineyd[nd, i] = (
                            (boundy(j + 1) - boundy(j)) * bound_locd[nd] / lenbound[j]
                        )
                    turbinex[i] = (
                        boundx(j)
                        + (boundx(j + 1) - boundx(j))
                        * (bound_loc - sum(lenbound[: j - 1]))
                        / lenbound[j]
                    )
                    turbiney[i] = (
                        boundy(j)
                        + (boundy(j + 1) - boundy(j))
                        * (bound_loc - sum(lenbound[: j - 1]))
                        / lenbound[j]
                    )
                    done = 1
                    bound_loc = bound_loc + spacing

    return turbinex, turbinexd, turbiney, turbineyd


def turbinelocations_dv(
    nbounds,
    nrows,
    nturbines,
    nouter,
    ngrid,
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
    start,
    startd,
    boundx,
    boundy,
    turbinex,
    turbinexd,
    turbiney,
    turbineyd,
    nbdirs,
):
    outerx = zeros(nouter)
    outery = zeros(nouter)
    outerxd = zeros((nbdirs, nouter))
    outeryd = zeros((nbdirs, nouter))
    innerx = zeros(ngrid)
    innery = zeros(ngrid)
    innerxd = zeros((nbdirs, ngrid))
    inneryd = zeros((nbdirs, ngrid))
    turbinex, turbinexd, turbiney, turbineyd = makeboundary_dv(
        nbounds,
        nouter,
        boundx,
        boundy,
        start,
        startd,
        outerx,
        outerxd,
        outery,
        outeryd,
        nbdirs,
    )
    turbinex, turbinexd, turbiney, turbineyd = makegrid_fortran_dv(
        nrows,
        ngrid,
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
        innerx,
        innerxd,
        innery,
        inneryd,
        nbdirs,
    )
    turbinexd[:, :] = 0.0_8
    turbineyd[:, :] = 0.0_8
    for nd in range(1, nbdirs):
        turbinexd[nd, 1:nouter] = outerxd[nd, 1]
        turbineyd[nd, 1:nouter] = outeryd[nd, 1]
        turbinexd[nd, nouter + 1 :] = innerxd[nd, 1]
        turbineyd[nd, nouter + 1 :] = inneryd[nd, 1]
    turbinex[1:nouter] = outerx[1]
    turbiney[1:nouter] = outery[1]
    turbinex[nouter + 1 :] = innerx[1]
    turbiney[nouter + 1 :] = innery[1]

    return turbinex, turbiney, turbinexd, turbineyd
