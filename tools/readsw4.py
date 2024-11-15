import numpy as np

def readimage(
    imfile: str,
    pnr: int, 
    verbose: bool = False
  ) -> tuple[np.ndarray, np.ndarray|np.float64, 
             np.ndarray|np.float64, np.ndarray|np.float64, int]:
  """ 
  Reads image produced by sw4 with *.sw4img extension.
  This script will not read the older WPP image format.

  im, x, y, z, readz = readimage(imfile, pnr, verbose)

  Parameters
  ----------
  imfile : string
      The path to the .sw4img file to be read.
  pnr : int
      The patch number to be read, e.g. for an image
      file spanning grids 0 to 2, to access the data
      for grid 1, pnr = 1.
  verbose : bool
      Set to true to display file header information.

  Returns
  -------
  im : np.ndarray
    The image, as a 2D array.
  x, y, z : np.ndarray | np.float64
    The spatial coordinates of the image. One of these is a scalar
    depending on which plane the image is in. The other two are the 
    same size as im.
  readz : int
    Grid/patch type indicator (0-cartesian, 1-curvilinear).
  """
  # Image mode string dictionary
  imgmodestr = ["none", "ux", "uy", "uz", "rho", "lambda",
                "mu", "p", "s", "uxexact", "uyexact", "uzexact",
                "div", "curlmag", "divdt", "curlmagdt", "lat",
                "lon", "topo", "gridx", "gridy", "gridz", "uxerr",
                "uyerr", "uzerr", "magdudt", "hmaxdudt", "vmaxdudt",
                "mag", "hmag", "hmax", "vmax", "gradrho", "gradmu",
                "gradlambda", "gradp", "grads", "qp", "qs"]

  # Open the file to be read
  with open(imfile, "rb") as fid:
    prec = np.fromfile(fid, dtype=np.int32, count=1)[0]
    npatches = np.fromfile(fid, dtype=np.int32, count=1)[0]
    t = np.fromfile(fid, dtype=np.float64, count=1)[0]
    plane = np.fromfile(fid, dtype=np.int32, count=1)[0]
    coord = np.fromfile(fid, dtype=np.float64, count=1)[0]

    mode = np.fromfile(fid, dtype=np.int32, count=1)[0]
    if mode < len(imgmodestr):
      mstr = imgmodestr[mode]
    else:
      mstr = 'unknown'

    gridinfo = np.fromfile(fid, dtype=np.int32, count=1)[0]
    timecreated = np.fromfile(fid, np.uint8, count=25)[:]

    # Display header
    if verbose:
      print(f"Found: prec = {prec}, t = {t}, plane = {plane}")
      print(f"       npatches = {npatches}, coord = {coord}, mode = {mstr}")
      print(f"       gridinfo (number of curvilinear patches) = {gridinfo}")
      print(f"       file created on {''.join(str(x) for x in timecreated)}")

    # Patch information
    firstCurviPatch = npatches - gridinfo
    h = np.zeros(npatches, dtype=np.float64)
    zmin = np.zeros(npatches, dtype=np.float64)
    ib = np.zeros(npatches, dtype=np.int32)
    ni = np.zeros(npatches, dtype=np.int32)
    jb = np.zeros(npatches, dtype=np.int32)
    nj = np.zeros(npatches, dtype=np.int32)
    for p in range(npatches):
      h[p] = np.fromfile(fid, dtype=np.float64, count=1)[0]
      zmin[p] = np.fromfile(fid, dtype=np.float64, count=1)[0]
      ib[p] = np.fromfile(fid, dtype=np.int32, count=1)[0]
      ni[p] = np.fromfile(fid, dtype=np.int32, count=1)[0]
      jb[p] = np.fromfile(fid, dtype=np.int32, count=1)[0]
      nj[p] = np.fromfile(fid, dtype=np.int32, count=1)[0]
      if verbose:
        print(f'Patch number {p} has h={h[p]}, zmin={zmin[p]:.2f}')

    ## Read data
    if pnr > (npatches-1):
      print(f'Error: input patch number (pnr={pnr}) must be in the range [0, {npatches-1}]')
    else:
      # Skip the preceding data patches
      for p in range(pnr):
        fid.seek( (ni[p]-ib[p]+1)*(nj[p]-jb[p]+1)*prec, 1 )
      # Read data for target patch pnr
      if prec == 4:
        im0 = np.fromfile(fid, dtype=np.float32, \
                          count=(ni[pnr]-ib[pnr]+1)*(nj[pnr]-jb[pnr]+1)) \
                          .reshape((ni[pnr]-ib[pnr]+1, nj[pnr]-jb[pnr]+1), order='F')
      else:
        im0 = np.fromfile(fid, dtype=np.float64, \
                          count=(ni[pnr]-ib[pnr]+1)*(nj[pnr]-jb[pnr]+1)) \
                          .reshape((ni[pnr]-ib[pnr]+1, nj[pnr]-jb[pnr]+1), order='F')
      # Skip the data patches after
      for p in range(pnr+1,npatches):
        fid.seek( (ni[p]-ib[p]+1)*(nj[p]-jb[p]+1)*prec, 1 )

      ## Read grid z-coordinates if curvilinear
      readz = 0
      # If target patch is a valid curvilinear patch
      if (pnr >= firstCurviPatch) and (gridinfo >= 1): 
        # Skip the preceding z-coordinate patches
        for p in range(firstCurviPatch, pnr):
          fid.seek( (ni[p]-ib[p]+1)*(nj[p]-jb[p]+1)*prec, 1 )
        # Read data for target z-coordinate patch
        if prec == 4:
          z0 = np.fromfile(fid, dtype=np.float32, \
                            count=(ni[pnr]-ib[pnr]+1)*(nj[pnr]-jb[pnr]+1)) \
                            .reshape((ni[pnr]-ib[pnr]+1, nj[pnr]-jb[pnr]+1), order='F')
        else:
          z0 = np.fromfile(fid, dtype=np.float64, \
                            count=(ni[pnr]-ib[pnr]+1)*(nj[pnr]-jb[pnr]+1)) \
                            .reshape((ni[pnr]-ib[pnr]+1, nj[pnr]-jb[pnr]+1), order='F')
        # Transpose z0 so it gets the same dimensions as im
        z = z0.T
        readz = 1

      ## Transpose im0 and return result in im
      im = im0.T

      ## Forming grid in x and y
      n1, n2 = im.shape
      if plane == 0: # y-z plane
        x = coord
        if readz == 0: # if pnr is cartesian patch
          # Make y and z two-dimensional (to match im in shape)
          y = np.zeros(im.shape)
          z = np.zeros(im.shape)
          for i in range(n1):
            y[i,:] = h[pnr]*(np.arange(ib[pnr], ni[pnr]+1)-1)
          for j in range(n2):
            z[:,j] = zmin[pnr] + h[pnr]*(np.arange(jb[pnr], nj[pnr]+1)-1)
        else: # curvilinear patch
          # Make y two-dimensional (to match z and im in shape)
          y = np.zeros(im.shape)
          for i in range(n1):
            y[i,:] = h[pnr]*(np.arange(ib[pnr], ni[pnr]+1)-1)
      elif plane == 1: # x-z plane
        y = coord
        if readz == 0: # if pnr is cartesian patch
          # Make x and z two-dimensional (to match im in shape)
          x = np.zeros(im.shape)
          z = np.zeros(im.shape)
          for i in range(n1):
            x[i,:] = h[pnr]*(np.arange(ib[pnr], ni[pnr]+1)-1)
          for j in range(n2):
            z[:,j] = zmin[pnr] + h[pnr]*(np.arange(jb[pnr], nj[pnr]+1)-1)
        else: # curvilinear patch
          # Make x two-dimensional (to match z and im in shape)
          x = np.zeros(im.shape)
          for i in range(n1):
            x[i,:] = h[pnr]*(np.arange(ib[pnr], ni[pnr]+1)-1)
      elif plane == 2: # x-y plane
        # Make x and y two-dimensional (to match im in shape)
        x = np.zeros(im.shape)
        y = np.zeros(im.shape)
        for i in range(n1):
          x[i,:] = h[pnr]*(np.arange(ib[pnr], ni[pnr]+1)-1)
        for j in range(n2):
          y[:,j] = h[pnr]*(np.arange(jb[pnr], nj[pnr]+1)-1)
        z = coord
      # end if pnr is valid patch
    # end open() as fid
  return im, x, y, z, readz

def readsac(
    filename : str, 
    verbose : bool = False
  ) -> tuple[np.ndarray,
             np.float32, np.float32, np.float32, np.float32, np.float32,
             np.float32, np.int32, np.int32, np.int32, np.int32, np.int32,
             np.int32, np.float32, np.float32, np.float32, str]:
  """ 
  Read SAC receiver data.

  u, dt, npts = readsac(filename, verbose)[0:3]

  u, dt, npts, t0, t1, lat, lon,
  year, jday, hour, min, sec, msec,
  cmpaz, cmpinc, idep, stnam = readsac(filename, verbose)

  Parameters
  ----------
  filename : string
      Name of the SAC file
  verbose : bool
      Set to true to display file header information.

  Returns
  -------
  u : np.ndarray
    The data component on the SAC file.
  dt : np.float32
    Uniform time step for u.
  npts : np.float32
    Number of elements in u.
  t0 : np.float32
    Begin time relative reference datum.
  t1 : np.float32
    End time relative reference datum.
  lat, lon : np.float32
    Latitude and longitude position of the SAC station.
  year, jday, hour, min, sec, msec : np.int32
    Reference data.
  cmpaz : np.float32
    Azimuth angle of component in degrees.
  cmpinc : np.float32
    Inclination angle of component in degrees.
  idep : np.int32
    Code of component stored (6-displacement, 7-velocity).
  stnam : str
    Name of the station.
  """
  # Open the file to be read
  with open(filename, "rb") as fid:
    dt = np.fromfile(fid, dtype=np.float32, count=1)[0]
    fid.seek( 4*4, 1 )
    t0 = np.fromfile(fid, dtype=np.float32, count=1)[0]
    t1 = np.fromfile(fid, dtype=np.float32, count=1)[0]
    fid.seek( 24*4, 1 )
    lat = np.fromfile(fid, dtype=np.float32, count=1)[0]
    lon = np.fromfile(fid, dtype=np.float32, count=1)[0]
    fid.seek( 2*4, 1 )
    evlat = np.fromfile(fid, dtype=np.float32, count=1)[0]
    evlon = np.fromfile(fid, dtype=np.float32, count=1)[0]
    fid.seek( 4, 1 )
    evdepth = np.fromfile(fid, dtype=np.float32, count=1)[0]
    fid.seek( 4*18, 1 )
    cmpaz  = np.fromfile(fid, dtype=np.float32, count=1)[0]
    cmpinc = np.fromfile(fid, dtype=np.float32, count=1)[0]
    fid.seek( 4*11, 1 )

    # Integers from offset 70
    year = np.fromfile(fid, dtype=np.int32, count=1)[0]
    jday = np.fromfile(fid, dtype=np.int32, count=1)[0]
    hour = np.fromfile(fid, dtype=np.int32, count=1)[0]
    min  = np.fromfile(fid, dtype=np.int32, count=1)[0]
    sec  = np.fromfile(fid, dtype=np.int32, count=1)[0]
    msec = np.fromfile(fid, dtype=np.int32, count=1)[0]
    nvhdr = np.fromfile(fid, dtype=np.int32, count=1)[0]
    fid.seek( 4*2, 1 )
    npts = np.fromfile(fid, dtype=np.int32, count=1)[0]
    fid.seek( 4*6, 1 )
    idep = np.fromfile(fid, dtype=np.int32, count=1)[0]
    fid.seek( 4*23, 1 )
    stnam = np.fromfile(fid, dtype=np.int8, count=8)
    stnam = ''.join([chr(item) for item in stnam])
    fid.seek( 4*46, 1 )

    if verbose:
      print(f'Station name = {stnam}')
      print(f'Year = {year}, Julian Day = {jday},')
      print(f'Hour = {hour}, Min = {min}, Sec = {sec}, Micro Sec = {msec},')
      print(f'Begin time (B) = {t0}, End time (E) = {t1}, dt = {dt}')
      print(f'Station (lat, lon) = ({lat}, {lon})')
      print(f'nvhdr = {nvhdr}, npts = {npts}')
      print(f'cmpaz = {cmpaz}, cmpinc = {cmpinc}, idep = {idep}')  

    # Read time series
    u = np.fromfile(fid, dtype=np.float32, count=npts)
    # end open() as fid
  return u, dt, npts, t0, t1, lat, lon, \
         year, jday, hour, min, sec, msec, cmpaz, cmpinc, idep, stnam