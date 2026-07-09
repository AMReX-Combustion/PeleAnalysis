.. highlight:: bash

makeImage
*********

Render 2-D images directly from AMReX plotfiles. For each plotfile and
variable, ``makeImage`` reads the finest-level data (the whole domain in 2-D,
or a single slice plane in 3-D), maps it through a colormap, and writes an
image file. It is a quick way to eyeball a field without opening VisIt or
ParaView, and to batch-generate frames for animations.

The tool is written to be **dependency-free**: it needs nothing beyond what the
rest of PeleAnalysis already links against. It can write three formats:

* **PPM** (``P6``, the default) — an uncompressed raw raster that every image
  viewer and converter understands.
* **PNG** — produced by a small self-contained encoder (no ``libpng``/``zlib``
  needed). The PNG pixel data is stored *uncompressed*, so the files are valid
  but larger than a typical PNG; see `Compressing the output`_ below to shrink
  them.
* **PDF** — a single page holding the raster, written by a small self-contained
  encoder. The image is embedded uncompressed and sized so one pixel maps to one
  PDF point; handy for dropping a field straight into a document.

Usage: ::

   ./makeImage2d.gnu.ex infile=PLT1 [PLT2 ...] vars="VAR1 VAR2 ..." [OPTIONS]

Example: ::

   ./makeImage2d.gnu.ex infile=plt00000 vars="temp Y(H2)" format=png
   ./makeImage2d.gnu.ex ./InputsSamples/makeImage.inp

Build: ::

   make EBASE=makeImage DIM=2       # 2-D
   make EBASE=makeImage DIM=3       # 3-D (renders one slice plane)

Tool Options
############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile    = plt00000 plt00100          # One or more plotfiles (space separated); required
   vars      = temp Y(H2)                 # Space-separated variable list; required
   outputDir = images                     # DEF: current dir; created if it does not exist
   format    = ppm                        # [ppm, png, pdf], DEF: ppm

``infile`` accepts a list of plotfiles; ``makeImage`` loops over all of them.
``vars`` is a space-separated list of variables. Both **stored** components and
**derived** quantities are accepted (anything ``amrData.CanDerive`` recognises);
an unknown name aborts the run. One image is written per plotfile/variable
combination, named ``<plotfile>_<var>[_<slice>].<ext>`` (any ``/`` in a derived
variable name is replaced by ``_``). For example ``plt00100`` and
``vars="temp"`` produce ``plt00100_temp.ppm``.
::

   #------------------- Colormap -------------------------------------------------------------
   colormap  = jet                        # DEF: jet; see the list below
   reverse   = 0                          # [0, 1], DEF: 0; reverse the colormap direction
   goPastMax = 1                          # [0, 1], DEF: 1; jet only (see below)

Two colormaps are built in: ``jet`` (the classic blue → cyan → green → yellow →
red ramp) and ``grayscale`` (low → black, high → white). In addition, a set of
Matplotlib colormaps is bundled (see `Colormaps`_ below). ``reverse=1`` flips
the direction of whichever colormap is selected.

``goPastMax`` only affects ``jet``. When ``1`` (default), data values *above*
the color-range maximum are not clipped to red but continue through magenta to
white, which makes over-range regions visually obvious. When ``0``, everything
above the maximum is drawn as a single dark red. ``goPastMax`` is ignored when
``reverse=1`` and for all colormaps other than ``jet``.
::

   #------------------- Color range ----------------------------------------------------------
   #useminmax1 = 300 2200                 # Fix [min max] for variable 1 (1-based index)
   #useminmax2 = 0 0.02                   # Fix [min max] for variable 2

By default each variable is normalised between its own data minimum and maximum
(computed over every level up to the level used). ``useminmax<i>`` fixes the
range for the *i*-th variable (1-based, in the order given in ``vars``), which
is what you want when generating a consistent series of frames across several
plotfiles. A constant field (min == max) is drawn at the bottom of the colormap
rather than dividing by zero.
::

   #------------------- Level / slice control ------------------------------------------------
   #finestLevel = 0                       # DEF: finest level in the file
   #yslice = 0                            # 3-D only; slice by cell index; DEF (3-D): yslice=0
   #xslice = 64                           # 3-D only; slice by cell index
   #zslice = 128                          # 3-D only; slice by cell index
   #ysliceCoord = 0.01                    # 3-D only; slice by physical coordinate (nearest cell)
   #xsliceCoord = -0.002                  # 3-D only; slice by physical coordinate (nearest cell)
   #zsliceCoord = -0.002                  # 3-D only; slice by physical coordinate (nearest cell)

``finestLevel`` caps the AMR level used to build the image (the data is sampled
on the uniform grid of that level); by default the finest level present in the
file is used, giving the highest-resolution picture.

In a **3-D** build the domain must be reduced to a plane. Specify exactly one
plane, either by **cell index** with ``xslice``/``yslice``/``zslice`` (index at
the level used) or by **physical coordinate** with
``xsliceCoord``/``ysliceCoord``/``zsliceCoord`` (the value is snapped to the
nearest cell of the level used); the default is ``yslice=0``. Giving more than
one plane aborts. Coordinates that fall outside the domain are clamped to the
nearest boundary cell. The two in-plane directions become the image width and
height, and the slice appears in the output filename, e.g.
``plt00100_temp_Y0.png``. In a **2-D** build these options are ignored and the
whole domain is rendered.

Colormaps
#########

Two colormaps are implemented directly in the tool:

* ``jet`` (default) — blue → cyan → green → yellow → red, with the optional
  ``goPastMax`` overshoot.
* ``grayscale`` — low → black, high → white.

In addition, a set of Matplotlib colormaps is bundled as lookup tables in
``Src/Colormaps.H`` (a small, dependency-free, header-only file that other
PeleAnalysis tools can reuse). Pass any of the names below to ``colormap=`` and
combine with ``reverse=1`` if you want the opposite direction:

* **Sequential:** ``Greys`` ``Purples`` ``Blues`` ``Greens`` ``Oranges``
  ``Reds`` ``YlOrBr`` ``YlOrRd`` ``OrRd`` ``PuRd`` ``RdPu`` ``BuPu`` ``GnBu``
  ``PuBu`` ``YlGnBu`` ``PuBuGn`` ``BuGn`` ``YlGn``
* **Sequential (2):** ``binary`` ``gist_yarg`` ``gist_gray`` ``gray`` ``bone``
  ``pink`` ``spring`` ``summer`` ``autumn`` ``winter`` ``cool`` ``Wistia``
  ``hot`` ``afmhot`` ``gist_heat`` ``copper``
* **Diverging:** ``PiYG`` ``PRGn`` ``BrBG`` ``PuOr`` ``RdGy`` ``RdBu``
  ``RdYlBu`` ``RdYlGn`` ``Spectral`` ``coolwarm`` ``bwr`` ``seismic``

These follow the Matplotlib groupings; see the `Matplotlib colormap reference
<https://matplotlib.org/stable/users/explain/colors/colormaps.html>`_ for what
each one looks like.

.. note::

   **Licensing.** The bundled color values were sampled from Matplotlib and are
   redistributed under their original licenses: the ColorBrewer maps (all
   *Sequential* maps plus ``PiYG``, ``PRGn``, ``BrBG``, ``PuOr``, ``RdGy``,
   ``RdBu``, ``RdYlBu``, ``RdYlGn`` and ``Spectral``) are © 2002 Cynthia Brewer,
   Mark Harrower and The Pennsylvania State University under the Apache License
   2.0; the remaining maps are © the Matplotlib Development Team under the
   Matplotlib license; ``coolwarm`` is by Kenneth Moreland. The full license
   texts and the map-to-license mapping are in ``licenses/``.

Image Orientation
#################

Images are written with the physical *+x* direction pointing right and *+y* (or
the second in-plane axis for a 3-D slice) pointing up, i.e. the data is flipped
vertically relative to raw array order so the picture is the right way up. Pixel
dimensions equal the number of cells of the level used, so a level-*L* image of
a domain with *N* base cells is *N* · 2\ :sup:`L` pixels across.

Output
######

One 8-bit RGB image per plotfile/variable, written to ``outputDir`` (default:
the current directory). PPM files can be opened directly by most image tools;
PNG and PDF files open anywhere. PNG and PDF embed the pixels uncompressed (see
below to shrink them).

Compressing the output
======================

The built-in PNG encoder stores pixels uncompressed to avoid pulling in a
compression library, so a ``.png`` from ``makeImage`` is roughly the size of the
equivalent raw image. If file size matters, run any standard optimiser as a
post-processing step — all of these rewrite the file as a normal compressed PNG:

.. code-block:: bash

   # Re-encode an existing PNG in place (lossless)
   optipng -o5 plt00100_temp.png
   pngcrush -ow plt00100_temp.png
   zopflipng -y plt00100_temp.png plt00100_temp.png

   # ImageMagick (also converts formats)
   magick plt00100_temp.png -strip plt00100_temp.png

If you output PPM instead, convert (and compress) it to PNG in one step:

.. code-block:: bash

   # Netpbm
   pnmtopng plt00100_temp.ppm > plt00100_temp.png

   # ImageMagick / GraphicsMagick
   magick plt00100_temp.ppm plt00100_temp.png
   convert   plt00100_temp.ppm plt00100_temp.png

   # Batch-convert a directory of frames
   for f in images/*.ppm; do magick "$f" "${f%.ppm}.png"; done

To assemble a series of frames into a movie once converted (for example with
``ffmpeg``):

.. code-block:: bash

   ffmpeg -framerate 24 -pattern_type glob -i 'images/plt*_temp.png' out.mp4

.. note::

   ``makeImage`` produces a single flat raster of the level used; it does not
   draw axes, colorbars, grid lines or AMR box outlines. For publication-quality
   figures with annotations, use VisIt or ParaView. ``makeImage`` is aimed at
   fast previews and bulk frame generation.
