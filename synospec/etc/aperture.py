"""
Define apertures to use for on-sky integrations.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import os
import numpy
import time

from scipy import signal

from shapely.geometry import Point, asPolygon
from shapely.affinity import rotate

def polygon_winding_number(polygon, point):
    """
    Determine the winding number of a 2D polygon about a point.  The
    code does **not** check if the polygon is simple (no interesecting
    line segments).  Algorithm taken from Numerical Recipies Section
    21.4.

    Args:
        polygon (`numpy.ndarray`_):
            An Nx2 array containing the x,y coordinates of a polygon.
            The points should be ordered either counter-clockwise or
            clockwise.
        point (`numpy.ndarray`_):
            One or more points for the winding number calculation.
            Must be either a 2-element array for a single (x,y) pair,
            or an Nx2 array with N (x,y) points.

    Returns:
        int or `numpy.ndarray`: The winding number of each point with
        respect to the provided polygon. Points inside the polygon
        have winding numbers of 1 or -1; see
        :func:`point_inside_polygon`.

    Raises:
        ValueError:
            Raised if `polygon` is not 2D, if `polygon` does not have
            two columns, or if the last axis of `point` does not have
            2 and only 2 elements.
    """
    # Check input shape is for 2D only
    if len(polygon.shape) != 2:
        raise ValueError('Polygon must be an Nx2 array.')
    if polygon.shape[1] != 2:
        raise ValueError('Polygon must be in two dimensions.')
    _point = numpy.atleast_2d(point)
    if _point.shape[1] != 2:
        raise ValueError('Point must contain two elements.')

    # Get the winding number
    nvert = polygon.shape[0]
    np = _point.shape[0]

    dl = numpy.roll(polygon, 1, axis=0)[None,:,:] - _point[:,None,:]
    dr = polygon[None,:,:] - point[:,None,:]
    dx = dl[:,:,0]*dr[:,:,1] - dl[:,:,1]*dr[:,:,0]

    indx_l = dl[:,:,1] > 0
    indx_r = dr[:,:,1] > 0

    wind = numpy.zeros((np, nvert), dtype=int)
    wind[indx_l & numpy.invert(indx_r) & (dx < 0)] = -1
    wind[numpy.invert(indx_l) & indx_r & (dx > 0)] = 1

    return numpy.sum(wind, axis=1)[0] if point.ndim == 1 else numpy.sum(wind, axis=1)


def point_inside_polygon(polygon, point):
    """
    Determine if one or more points is inside the provided polygon.

    Primarily a wrapper for :func:`polygon_winding_number`, that
    returns True for each poing that is inside the polygon.

    Args:
        polygon (`numpy.ndarray`_):
            An Nx2 array containing the x,y coordinates of a polygon.
            The points should be ordered either counter-clockwise or
            clockwise.
        point (`numpy.ndarray`_):
            One or more points for the winding number calculation.
            Must be either a 2-element array for a single (x,y) pair,
            or an Nx2 array with N (x,y) points.

    Returns:
        bool or `numpy.ndarray`: Boolean indicating whether or not
        each point is within the polygon.
    """
    return numpy.absolute(polygon_winding_number(polygon, point)) == 1


def hexagon_vertices(d=1., incircle=False):
    r"""
    Construct the vertices of a hexagon.

    The long axis of the hexagon is always oriented along the Cartesian
    :math:`x` axis.

    Args:
        d (:obj:`float`, optional):
            The diameter of circumscribed circle.
        incircle (:obj:`bool`, optional):
            Use the provided diameter to set the incircle of the hexagon
            instead of its circumscribed circle.

    Returns:
        `numpy.ndarray`_: An array with shape :math:`(6,2)`, providing the x
        and y Cartesian vertices of the hexagon.
    """
    # Get the incircle radius
    r = d/2/numpy.cos(numpy.pi/6) if incircle else d/2 
    # Generate each vertex, in clockwise order, using a brute force approach.
    v = numpy.zeros((6,2), dtype=float)
    sina = numpy.sin(numpy.pi/3.0)

    v[0,0] = -r/2
    v[1,0] = r/2
    v[2,0] = r
    v[3,0] = r/2
    v[4,0] = -r/2
    v[5,0] = -r

    v[0,1] = r * sina
    v[1,1] = r * sina
    v[2,1] = 0.
    v[3,1] = -r * sina
    v[4,1] = -r * sina
    v[5,1] = 0.
    
    return v


class Aperture:
    """
    Abstract class for a general aperture shape.

    .. todo:
        - limit the calculation of the polygon overlap to where the
        corners of the grid cells cross the boundary of the aperture.

    Args:
        shape (`shapely.geometry.base.BaseGeometry`_):
            A shape object from the Shapely python package.
    """
    def __init__(self, shape):
#        self.shape = shapely.prepared.prep(shape)
        self.shape = shape

    def response(self, x, y, method='fractional'):
        r"""
        Compute the response function of the aperture to the sky over
        a regular grid. This is the same as rendering an "image" of
        the aperture.

        The integral of the returned map is normalized to the
        aperture area.

        Args:
            x (array-like):
                The list of x coordinates for the grid.  Must be
                linearly spaced.
            y (array-like):
                The list of y coordinates for the grid.  Must be
                linearly spaced.
            method (:obj:`str`, optional):
                Method used to construct the overlap grid.  Options
                are:

                    - 'whole': Any grid cell with its center inside the
                      aperture is set to the area of the grid cell.  All
                      others set to 0.
                    - 'fractional': Perform the detailed calculation of
                      the fraction of each grid-cell within the
                      aperture.
    
        Returns:
            `numpy.ndarray`_: An array with shape :math:`(N_x, N_y)`
            with the fraction of each grid cell covered by the
            aperture.

        Raises:
            ValueError:
                Raised if the provided arguments are not regularly
                spaced, or if there aren't at least 2 grid points in
                each dimension.
        """
        # Check input
        if len(x) < 2 or len(y) < 2:
            raise ValueError('Must provide at least 2 points per grid point.')
        minimum_x_difference = 0 if numpy.issubdtype(x.dtype, numpy.integer) \
                                    else numpy.finfo(x.dtype).eps*100
        minimum_y_difference = 0 if numpy.issubdtype(y.dtype, numpy.integer) \
                                    else numpy.finfo(y.dtype).eps*100
        if numpy.any(numpy.absolute(numpy.diff(numpy.diff(x))) > minimum_x_difference):
            raise ValueError('X coordinates are not regular to numerical precision.')
        if numpy.any(numpy.absolute(numpy.diff(numpy.diff(y))) > minimum_y_difference):
            raise ValueError('Y coordinates are not regular to numerical precision.')

        # Grid shape
        nx = len(x)
        ny = len(y)

        # Get the cell size
        dx = abs(x[1]-x[0])
        dy = abs(y[1]-y[0])

        cell_area = dx*dy

        if method == 'whole':
            # Only include whole pixels
#            X,Y = map(lambda x : x.ravel(), numpy.meshgrid(x, y))
#            img = numpy.array(list(map(lambda x: self.shape.contains(Point(x[0],x[1])),
#                                        zip(X,Y)))).reshape(ny,nx).astype(int)/cell_area
            # Build the coordinate of the pixel centers
            coo = numpy.array(list(map(lambda x : x.ravel(), numpy.meshgrid(x, y)))).T
            # Find the pixels with their centers inside the shape
            img = point_inside_polygon(self.vertices(), coo).reshape(ny,nx).astype(int)/cell_area
            # Return after adjusting to match the defined area
            return img * (self.area/numpy.sum(img)/cell_area)
        elif method == 'fractional':
            # Allow for fractional pixels by determining the overlap
            # between the shape and each grid cell

#            # Build the cell polygons
#            cells, sx, ex, sy, ey = self._overlapping_grid_polygons(x, y)
#
#            # Construct a grid with the fractional area covered by the
#            # aperture
#            img = numpy.zeros((len(y), len(x)), dtype=float)
#            img[sy:ey,sx:ex] = numpy.array(list(map(lambda x: self.shape.intersection(x).area,
#                                                    cells))).reshape(ey-sy,ex-sx)/cell_area
#            return img * (self.area/numpy.sum(img)/cell_area)

            # NOTE: This is much faster than the above
            sx, ex, _, sy, ey, _ = self._overlapping_region(x,y)

            X,Y = map(lambda x : x.ravel(), numpy.meshgrid(x[sx:ex], y[sy:ey]))

            # Get the points describing the corners of each overlapping grid cell
            cx = X[:,None] + (numpy.array([-0.5,0.5,0.5,-0.5])*dx)[None,:]
            cy = Y[:,None] + (numpy.array([-0.5,-0.5,0.5,0.5])*dy)[None,:]
            cells = numpy.append(cx, cy, axis=1).reshape(-1,2,4).transpose(0,2,1).reshape(-1,2)
            # cells has shape (ncell*4,2)
            inside_shape = point_inside_polygon(self.vertices(), cells).reshape(-1,4)
            indx = numpy.all(inside_shape, axis=1)
            _img = indx.astype(float)
            intersect = numpy.any(inside_shape, axis=1) & numpy.invert(indx)
            _img[intersect] = numpy.array(list(map(lambda x:
                                                    self.shape.intersection(asPolygon(x)).area,
                                                cells.reshape(-1,4,2)[intersect,...])))/cell_area
            img = numpy.zeros((len(y), len(x)), dtype=float)
            img[sy:ey,sx:ex] = _img.reshape(ey-sy,ex-sx)
            return img * (self.area/numpy.sum(img)/cell_area)

        raise ValueError('Unknown response method {0}.'.format(method))

#        # OLD and slow
#        cells, tree = Aperture._get_grid_tree(x, y, fast=False)
#        alpha = numpy.zeros((nx,ny), dtype=float)
#        for k in tree.intersection(self.shape.bounds):
#            i = k//ny
#            j = k - i*ny
#            alpha[i,j] = cells[k].intersection(self.shape).area
#        return alpha

    def _overlapping_region(self, x, y):
        r"""
        Return the starting and ending indices of the grid cells that
        overlap with the shape bounds.

        Args:
            x (array-like):
                The list of x coordinates for the grid.  Must be
                linearly spaced.
            y (array-like):
                The list of y coordinates for the grid.  Must be
                linearly spaced.
        
        Returns:
            :obj:`tuple`: Returns six scalars: the starting and
            ending x index, the grid step in x, the starting and
            ending y index, and the grid step in y for the region of
            the grid that overlaps the shape boundary.
        """
        # Get the cell size
        dx = abs(x[1]-x[0])
        dy = abs(y[1]-y[0])

        # Find the x coordinates of the grid cells that overlap the shape
        xlim = list(self.shape.bounds[::2])
        if xlim[0] > xlim[1]:
            xlim = xlim[::-1]
        xindx = (x+0.5*dx > xlim[0]) & (x-0.5*dx < xlim[1])
        sx = numpy.arange(len(x))[xindx][0]
        ex = numpy.arange(len(x))[xindx][-1]+1

        # Find the y coordinates of the grid cells that overlap the shape
        ylim = list(self.shape.bounds[1::2])
        if ylim[0] > ylim[1]:
            ylim = ylim[::-1]
        yindx = (y+0.5*dy > ylim[0]) & (y-0.5*dy < ylim[1])
        sy = numpy.arange(len(y))[yindx][0]
        ey = numpy.arange(len(y))[yindx][-1]+3

        return sx, ex, dx, sy, ey, dy

    def _overlapping_grid_polygons(self, x, y):
        r"""
        Construct a list of grid-cell polygons (rectangles) that are
        expected to overlap the aperture.

        The list of polygons follows array index order.  I.e., polygon
        :math:`k` is the cell at location :math:`(j,i)`, where::

        .. math::
            
            j = k//nx
            i = k - j*nx

        Args:
            x (array-like):
                The list of x coordinates for the grid.  Must be
                linearly spaced.
            y (array-like):
                The list of y coordinates for the grid.  Must be
                linearly spaced.
        
        Returns:
            :obj:`tuple`: Five objects are returned:

                - A list of `shapely.geometry.polygon.Polygon`_
                  objects, one per grid cell. Only those grid cells
                  that are expected to overlap the shape's bounding
                  box are included.
                - The starting and ending x index and the starting
                  and ending y index for the returned list of cell
                  polygons.

        """
        sx, ex, dx, sy, ey, dy = self._overlapping_region(x, y)

        # Construct the grid
        X,Y = map(lambda x : x.ravel(), numpy.meshgrid(x[sx:ex], y[sy:ey]))

        # Construct the polygons
        cx = X[:,None] + (numpy.array([-0.5,0.5,0.5,-0.5])*dx)[None,:]
        cy = Y[:,None] + (numpy.array([-0.5,-0.5,0.5,0.5])*dy)[None,:]

        boxes = numpy.append(cx, cy, axis=1).reshape(-1,2,4).transpose(0,2,1)
        polygons = [asPolygon(box) for box in boxes]
        
        return polygons, sx, ex, sy, ey

    @property
    def area(self):
        """The area of the aperture."""
        return self.shape.area

    @property
    def bounds(self):
        """The bounding box of the aperture."""
        return self.shape.bounds

    def integrate_over_source(self, source, response_method='fractional', sampling=None,
                              size=None):
        """
        Integrate the flux of a source over the aperture.

        This is done by generating an image of the aperture over the map
        of the source surface-brightness distribution, using
        :func:`Aperture.response`.  The source is expected to already
        have been mapped using its `make_map` function, or one should
        provide `sampling` and `size` values to construct the map inside
        this function.

        See also: :func:`Aperture.map_integral_over_source`.

        .. todo::

            No type checking is done to require that ``source`` is a
            :class:`~synospec.etc.source.Source` object, but the code
            will barf if it isn't.

        Args:
            source (:class:`~synospec.etc.source.Source`):
                Source surface-brightness distribution
            response_method (:obj:`str`, optional):
                See ``method`` argument for
                :func:`Aperture.response`.
            sampling (:obj:`float`, optional):
                Sampling of the square map in arcsec/pixel. If not
                None, the source map is (re)constructed.
            size (:obj:`float`, optional):
                Size of the square map in arcsec.  If not None, the
                source map is (re)constructed.

        Returns:
            :obj:`float`: The integral of the source over the
            aperture.

        Raises:
            ValueError:
                Raised if the source map has not been constructed and
                ``sampling`` and ``size`` are both None.
        """
        if source.data is None and sampling is None and size is None:
            raise ValueError('Must make a map of the source first.')

        if sampling is not None or size is not None:
            source.make_map(sampling=sampling, size=size)

        aperture_image = self.response(source.x, source.y, method=response_method)
        return numpy.sum(source.data*aperture_image)*numpy.square(source.sampling)

    def map_integral_over_source(self, source, response_method='fractional', sampling=None,
                                 size=None):
        """
        Construct a continuous map of the source integrated over the
        aperture.

        This is done by generating an image of the aperture over the map
        of the source surface-brightness distribution, using
        :func:`Aperture.response`.  The integral of the source over the
        aperture *at any offset position within the map* is calculated
        by convolving the the source distribution and the aperture
        image.

        See also :func:`Aperture.integrate_over_source`. A single
        call to this function or
        :func:`Aperture.integrate_over_source` to get the integral
        with no offset of the aperture are marginally different.
        However, use of this function is much more efficient if you
        want to calculate the integral of the source over many
        positional offsets of the aperture.
        
        .. todo::

            No type checking is done to require that ``source`` is a
            :class:`~synospec.etc.source.Source` object, but the code
            will barf if it isn't.

        Args:

            source (:class:`~synospec.etc.source.Source`):
                Source surface-brightness distribution
            response_method (:obj:`str`, optional):
                See ``method`` argument for
                :func:`Aperture.response`.
            sampling (:obj:`float`, optional):
                Sampling of the square map in arcsec/pixel. If not
                None, the source map is (re)constructed.
            size (:obj:`float`, optional):
                Size of the square map in arcsec.  If not None, the
                source map is (re)constructed.

        Returns:
            `numpy.ndarray`_: The integral of the source over the
            aperture with the aperture centered at any position in
            the map. The integral with no offset between the image of
            the aperture and the image of the source is::

                cy = source.data.shape[0]//2
                cx = source.data.shape[1]//2
                integral = Aperture.map_integral_over_source(source)[cy,cx]

            which should be identical to::

                integral = Aperture.integrate_over_source(source)

        """
        if source.data is None and sampling is None and size is None:
            raise ValueError('Must make a map of the source first.')

        if sampling is not None or size is not None:
            source.make_map(sampling=sampling, size=size)

        aperture_image = self.response(source.x, source.y, method=response_method)
        return signal.fftconvolve(source.data, aperture_image*numpy.square(source.sampling),
                                  mode='same')

    def vertices(self, wrap=False):
        vert = numpy.array(self.shape.exterior.coords)
        return vert if wrap else vert[:-1]
        

class FiberAperture(Aperture):
    """
    Define a fiber aperture.

    Note that the units for the center and diameter are only relevant
    in the application of the aperture to a source. They should
    typically be in arcseconds, with the center being relative to the
    source to observe.

    Args:
        cx (scalar-like):
            Center X coordinate, typically 0.
        cy (scalar-like):
            Center Y coordinate, typically 0.
        d (scalar-like):
            Fiber diameter.  Aperture is assumed to be a circle resolved
            by a set of line segments.
        resolution (:obj:`int`, optional):
            Set the "resolution" of the circle. Higher numbers mean
            more line segments are used to define the circle, but
            there isn't a 1-1 correspondence. See `shapely.buffer`_.
            Default is to use the `shapely`_ default.

    Attributes:
        center (:obj:`list`):
            Center x and y coordinate.
        diameter (:obj:`float`):
            Fiber diameter
    """
    def __init__(self, cx, cy, d, resolution=None):
        self.center = [cx,cy]
        self.diameter = d
        kw = {} if resolution is None else {'resolution':resolution}
        super(FiberAperture, self).__init__(Point(cx,cy).buffer(d/2, **kw))


# A hack the just creates a pseudonym for the FiberAperture
class CircularAperture(FiberAperture):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class SlitAperture(Aperture):
    """
    Define a slit aperture.

    The orientation of the slit is expected to have the length along
    the y axis and the width along the x axis. The rotation is
    counter-clockwise in a right-handed Cartesian frame.

    Note that the units for the center, width, and length are only
    relevant in the application of the aperture to a source. They
    should typically be in arcseconds, with the center being relative
    to the source to observe.

    Exactly the same aperture is obtained in the following two
    calls::

        s = SlitAperture(0., 0., 1, 10)
        ss = SlitAperture(0., 0., 10, 1, rotation=90)

    Args:
        cx (scalar-like):
            Center X coordinate, typically 0.
        cy (scalar-like):
            Center Y coordinate, typically 0.
        width (scalar-like):
            Slit width along the unrotated x axis.
        length (scalar-like):
            Slit length along the unrotated y axis.
        rotation (scalar-like):
            Cartesian rotation of the slit in degrees.

    Attributes:
        center (:obj:`list`):
            Center x and y coordinate.
        width (:obj:`float`):
            Slit width
        length (:obj:`float`):
            Slit length
        rotation (:obj:`float`):
            Slit rotation (deg)
    """
    def __init__(self, cx, cy, width, length, rotation=0.):
        self.center = [cx,cy]
        self.width = width
        self.length = length
        self.rotation = rotation
        x = numpy.array([-width/2, width/2])+cx
        y = numpy.array([-length/2, length/2])+cy
        square = asPolygon(numpy.append(numpy.roll(numpy.repeat(x,2),-1),
                                        numpy.repeat(y,2)).reshape(2,4).T)
        # rotate() function is provided by shapely.affinity package
        super(SlitAperture, self).__init__(rotate(square, rotation))


class HexagonalAperture(Aperture):
    """
    Define a hexagonal aperture.

    Note that the units for the center and diameter are only relevant in the
    application of the aperture to a source. They should typically be in
    arcseconds, with the center being relative to the source to observe.

    Args:
        cx (scalar-like):
            Center X coordinate, typically 0.
        cy (scalar-like):
            Center Y coordinate, typically 0.
        d (:obj:`float`):
            The diameter of circumscribed circle.
        incircle (:obj:`bool`, optional):
            Use the provided diameter to set the incircle of the hexagon
            instead of its circumscribed circle.
        orientation (:obj:`str`, :obj:`float`, optional):
            Sets the orientation of the hexagon, must be either
            'horizontal', 'vertical', or a rotation angle in degrees
            relative to the horizontal orientation.  The horizontal and
            vertical orientations set the long axis of the hexagon along
            the Cartesian x and y axis, respectively.  The horizontal
            orientation is equivalent to a rotation angle of 0 and a
            vertical orientation is equivalent to a rotation angle of 30
            degrees.  While the polar-coordinate ordering of the
            vertices in the output array will change, note the shape
            is degenerate with a periodicity of 60 degrees.

    Attributes:
        center (:obj:`list`):
            Center x and y coordinate.
        width (:obj:`float`):
            Slit width
        length (:obj:`float`):
            Slit length
        rotation (:obj:`float`):
            Slit rotation (deg)
    """
    def __init__(self, cx, cy, d, incircle=False, orientation='horizontal'):

        self.center = numpy.array([cx,cy])
        if orientation == 'horizontal':
            self.rotation = 0.
        elif orientation == 'vertical':
            self.rotation = 90.
        elif not isinstance(orientation, (int, float, numpy.integer, numpy.floating)):
            self.rotation = orientation
        else:
            raise ValueError('Orientation must be "horizontal", "vertical", or a numerical '
                             f'rotation in degrees.  Cannot interpret {orientation}.')

        hexv = hexagon_vertices(d=d, incircle=incircle)
        hexv[:,0] += cx
        hexv[:,1] += cy
        hexv = asPolygon(hexv)
        # rotate() function is provided by shapely.affinity package
        super().__init__(rotate(hexv, self.rotation))


