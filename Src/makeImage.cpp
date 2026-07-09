// makeImage: render 2-D images from AMReX plotfiles.
//
// For each requested plotfile / variable, the finest-level data (a 2-D domain,
// or a single slice plane of a 3-D domain) is mapped through a colormap and
// written as an image. Two output formats are supported:
//
//   * PPM (P6, the default) -- a trivial raw raster that any image viewer or
//     converter understands (see the documentation for one-line conversions
//     to PNG).
//   * PNG -- written by a small, self-contained encoder (no external image
//     library is required). The PNG data is stored uncompressed; the
//     documentation shows how to shrink the files with a standard optimizer.
//
// The tool is deliberately dependency-free so that it builds anywhere the rest
// of PeleAnalysis builds.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include <AMReX_DataServices.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <Colormaps.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Usage:\n"
    << "  " << argv[0]
    << " infile=PLT1 [PLT2 ...] vars=\"VAR1 VAR2 ...\" "
       "[OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=PLT1 [PLT2 ...]  One or more AMReX plotfiles to read\n"
    << "  vars=\"VAR1 ...\"         Space-separated list of variables to "
       "render\n\n"

    << "Options:\n"
    << "  format=ppm|png|pdf      Output image format (DEF: ppm)\n"
    << "  colormap=NAME           jet, grayscale, or a matplotlib colormap "
       "(see docs); DEF: jet\n"
    << "  reverse=0|1             Reverse the colormap direction (DEF: 0)\n"
    << "  goPastMax=0|1           jet only: extend past the max into "
       "magenta/white (DEF: 1)\n"
    << "  useminmax<i>=MIN MAX    Fix the color range for variable i "
       "(1-based)\n"
    << "  finestLevel=N           Cap the AMR level used (DEF: finest in "
       "file)\n"
    << "  outputDir=DIR           Directory for the images (DEF: current "
       "dir)\n"
#if (AMREX_SPACEDIM == 3)
    << "  xslice=I|yslice=J|zslice=K  3-D slice plane by cell index at finest "
       "level (DEF: yslice=0)\n"
    << "  xsliceCoord=X|ysliceCoord=Y|zsliceCoord=Z  3-D slice plane by "
       "physical coordinate (nearest cell)\n"
#endif
    << "  verbose                 Enable verbose plotfile I/O\n"
    << "  -h, --help              Show this help message\n\n"

    << "Output files are named <plotfile>_<var>[_<slice>].<ext>.\n";

  std::exit(1);
}

// Replace characters that are awkward in filenames (e.g. the '/' in derived
// names such as "x_velocity/density") with underscores.
static std::string
ProtectSlashes(const std::string& str)
{
  std::string s = str;
  for (char& c : s)
    if (c == '/')
      c = '_';
  return s;
}

// Basename of a (possibly trailing-slash-terminated) path, e.g.
// "../run/plt00100/" -> "plt00100".
static std::string
PlotfileBasename(const std::string& path)
{
  std::string s = path;
  while (s.size() > 1 && s.back() == '/')
    s.pop_back();
  const size_t pos = s.find_last_of('/');
  return (pos == std::string::npos) ? s : s.substr(pos + 1);
}

static unsigned char
ToByte(Real x)
{
  int v = (int)(x + 0.5);
  if (v < 0)
    v = 0;
  if (v > 255)
    v = 255;
  return (unsigned char)v;
}

// The classic PeleAnalysis "jet" ramp: blue -> cyan -> green -> yellow -> red,
// optionally continuing past the maximum (c > 1) through magenta to white.
static void
ColormapJet(
  Real c, bool goPastMax, unsigned char& r, unsigned char& g, unsigned char& b)
{
  if (c < 0.125) {
    r = 0;
    g = 0;
    b = ToByte((c + 0.125) * 1020.0);
  } else if (c < 0.375) {
    r = 0;
    g = ToByte((c - 0.125) * 1020.0);
    b = 255;
  } else if (c < 0.625) {
    r = ToByte((c - 0.375) * 1020.0);
    g = 255;
    b = ToByte((0.625 - c) * 1020.0);
  } else if (c < 0.875) {
    r = 255;
    g = ToByte((0.875 - c) * 1020.0);
    b = 0;
  } else if (c < 1.000) {
    r = ToByte((1.125 - c) * 1020.0);
    g = 0;
    b = 0;
  } else if (goPastMax) {
    if (c < 1.125) {
      r = ToByte((c - 0.875) * 1020.0);
      g = 0;
      b = ToByte((c - 1.000) * 1020.0);
    } else if (c < 1.250) {
      r = 255;
      g = 0;
      b = ToByte((c - 1.000) * 1020.0);
    } else if (c < 1.500) {
      r = 255;
      g = ToByte((c - 1.250) * 1020.0);
      b = 255;
    } else {
      r = 255;
      g = 255;
      b = 255;
    }
  } else {
    r = 128;
    g = 0;
    b = 0;
  }
}

static void
ColormapGray(Real c, unsigned char& r, unsigned char& g, unsigned char& b)
{
  r = g = b = ToByte(c * 255.0);
}

// -------------------------------------------------------------------------
// Image writers. Both take an 8-bit RGB buffer of width*height*3 bytes, row
// major, with the top row first.
// -------------------------------------------------------------------------

static bool
WritePPM(
  const std::string& path,
  int width,
  int height,
  const std::vector<unsigned char>& rgb)
{
  FILE* fp = fopen(path.c_str(), "wb");
  if (fp == nullptr)
    return false;
  fprintf(fp, "P6\n%d %d\n255\n", width, height);
  fwrite(rgb.data(), 1, (size_t)width * height * 3, fp);
  fclose(fp);
  return true;
}

// --- Minimal, dependency-free PNG support (8-bit RGB) ---

static uint32_t
Crc32(const unsigned char* buf, size_t len)
{
  static uint32_t table[256];
  static bool have_table = false;
  if (!have_table) {
    for (uint32_t n = 0; n < 256; ++n) {
      uint32_t c = n;
      for (int k = 0; k < 8; ++k)
        c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
      table[n] = c;
    }
    have_table = true;
  }
  uint32_t crc = 0xFFFFFFFFu;
  for (size_t i = 0; i < len; ++i)
    crc = table[(crc ^ buf[i]) & 0xFFu] ^ (crc >> 8);
  return crc ^ 0xFFFFFFFFu;
}

static uint32_t
Adler32(const unsigned char* buf, size_t len)
{
  const uint32_t MOD = 65521u;
  uint32_t a = 1, b = 0;
  size_t i = 0;
  while (i < len) {
    const size_t chunk = (len - i < 5552) ? (len - i) : 5552;
    for (size_t j = 0; j < chunk; ++j) {
      a += buf[i + j];
      b += a;
    }
    a %= MOD;
    b %= MOD;
    i += chunk;
  }
  return (b << 16) | a;
}

static void
PushU32BE(std::vector<unsigned char>& v, uint32_t x)
{
  v.push_back((unsigned char)((x >> 24) & 0xFF));
  v.push_back((unsigned char)((x >> 16) & 0xFF));
  v.push_back((unsigned char)((x >> 8) & 0xFF));
  v.push_back((unsigned char)(x & 0xFF));
}

static void
WritePNGChunk(
  FILE* fp, const char* type, const std::vector<unsigned char>& data)
{
  std::vector<unsigned char> len_be;
  PushU32BE(len_be, (uint32_t)data.size());
  fwrite(len_be.data(), 1, 4, fp);
  fwrite(type, 1, 4, fp);
  if (!data.empty())
    fwrite(data.data(), 1, data.size(), fp);
  // CRC covers the chunk type followed by the chunk data.
  std::vector<unsigned char> crc_in;
  crc_in.insert(crc_in.end(), type, type + 4);
  crc_in.insert(crc_in.end(), data.begin(), data.end());
  std::vector<unsigned char> crc_be;
  PushU32BE(crc_be, Crc32(crc_in.data(), crc_in.size()));
  fwrite(crc_be.data(), 1, 4, fp);
}

// Write an 8-bit RGB PNG using only uncompressed (stored) DEFLATE blocks, so
// no compression library is needed. The resulting file is a valid PNG; it is
// simply larger than a compressed one (see the docs for optimization tips).
static bool
WritePNG(
  const std::string& path,
  int width,
  int height,
  const std::vector<unsigned char>& rgb)
{
  FILE* fp = fopen(path.c_str(), "wb");
  if (fp == nullptr)
    return false;

  const unsigned char sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
  fwrite(sig, 1, 8, fp);

  std::vector<unsigned char> ihdr;
  PushU32BE(ihdr, (uint32_t)width);
  PushU32BE(ihdr, (uint32_t)height);
  ihdr.push_back(8); // bit depth
  ihdr.push_back(2); // color type: truecolour (RGB)
  ihdr.push_back(0); // compression method
  ihdr.push_back(0); // filter method
  ihdr.push_back(0); // interlace method
  WritePNGChunk(fp, "IHDR", ihdr);

  // Filtered image data: each scanline is prefixed with a filter-type byte (0).
  const size_t row_bytes = (size_t)width * 3;
  std::vector<unsigned char> raw;
  raw.reserve((row_bytes + 1) * (size_t)height);
  for (int y = 0; y < height; ++y) {
    raw.push_back(0); // filter: None
    raw.insert(
      raw.end(), rgb.begin() + (size_t)y * row_bytes,
      rgb.begin() + (size_t)y * row_bytes + row_bytes);
  }

  // zlib stream: 2-byte header + stored DEFLATE blocks + 4-byte Adler32.
  std::vector<unsigned char> zlib;
  zlib.push_back(0x78); // CMF
  zlib.push_back(0x01); // FLG (checksum ok, no dict, level 0)
  const size_t MAXBLK = 65535;
  const size_t total = raw.size();
  size_t off = 0;
  if (total == 0) {
    const unsigned char empty_block[5] = {0x01, 0x00, 0x00, 0xFF, 0xFF};
    zlib.insert(zlib.end(), empty_block, empty_block + 5);
  }
  while (off < total) {
    const size_t blk = (total - off < MAXBLK) ? (total - off) : MAXBLK;
    const unsigned char bfinal = (off + blk >= total) ? 1 : 0;
    const unsigned char len_lo = (unsigned char)(blk & 0xFF);
    const unsigned char len_hi = (unsigned char)((blk >> 8) & 0xFF);
    zlib.push_back(bfinal); // BFINAL in bit 0, BTYPE=00 (stored)
    zlib.push_back(len_lo);
    zlib.push_back(len_hi);
    zlib.push_back((unsigned char)~len_lo);
    zlib.push_back((unsigned char)~len_hi);
    zlib.insert(zlib.end(), raw.begin() + off, raw.begin() + off + blk);
    off += blk;
  }
  PushU32BE(zlib, Adler32(raw.data(), raw.size()));
  WritePNGChunk(fp, "IDAT", zlib);

  const std::vector<unsigned char> empty;
  WritePNGChunk(fp, "IEND", empty);

  fclose(fp);
  return true;
}

// Write a minimal single-page PDF that embeds the image as an uncompressed
// 8-bit RGB image XObject (top row first). No external library is needed; the
// page is sized so one pixel maps to one PDF point.
static bool
WritePDF(
  const std::string& path,
  int width,
  int height,
  const std::vector<unsigned char>& rgb)
{
  std::vector<unsigned char> buf;
  auto put = [&](const std::string& s) {
    buf.insert(buf.end(), s.begin(), s.end());
  };
  auto pad10 = [](size_t v) {
    std::string s = std::to_string(v);
    return std::string(s.size() < 10 ? 10 - s.size() : 0, '0') + s;
  };

  put("%PDF-1.4\n");
  put("%\xE2\xE3\xCF\xD3\n"); // binary-content marker

  const int nobj = 5;
  std::vector<size_t> off;

  off.push_back(buf.size());
  put("1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n");

  off.push_back(buf.size());
  put("2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n");

  off.push_back(buf.size());
  put(
    "3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 " +
    std::to_string(width) + " " + std::to_string(height) +
    "] /Resources << /XObject << /Im0 4 0 R >> >> /Contents 5 0 R >>\n"
    "endobj\n");

  off.push_back(buf.size());
  {
    const size_t n = (size_t)width * height * 3;
    put(
      "4 0 obj\n<< /Type /XObject /Subtype /Image /Width " +
      std::to_string(width) + " /Height " + std::to_string(height) +
      " /ColorSpace /DeviceRGB /BitsPerComponent 8 /Length " +
      std::to_string(n) + " >>\nstream\n");
    buf.insert(buf.end(), rgb.begin(), rgb.end());
    put("\nendstream\nendobj\n");
  }

  off.push_back(buf.size());
  {
    const std::string content = "q\n" + std::to_string(width) + " 0 0 " +
                                std::to_string(height) +
                                " 0 0 cm\n/Im0 Do\nQ\n";
    put(
      "5 0 obj\n<< /Length " + std::to_string(content.size()) +
      " >>\nstream\n" + content + "endstream\nendobj\n");
  }

  const size_t xrefPos = buf.size();
  put("xref\n0 " + std::to_string(nobj + 1) + "\n");
  put("0000000000 65535 f \n");
  for (int i = 0; i < nobj; ++i)
    put(pad10(off[i]) + " 00000 n \n");
  put(
    "trailer\n<< /Size " + std::to_string(nobj + 1) +
    " /Root 1 0 R >>\nstartxref\n" + std::to_string(xrefPos) + "\n%%EOF\n");

  FILE* fp = fopen(path.c_str(), "wb");
  if (fp == nullptr)
    return false;
  fwrite(buf.data(), 1, buf.size(), fp);
  fclose(fp);
  return true;
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    if (
      argc < 2 || std::string(argv[1]) == "-h" ||
      std::string(argv[1]) == "--help")
      print_usage(argc, argv);

    ParmParse pp;
    if (pp.contains("help"))
      print_usage(argc, argv);

    const bool verbose = ParallelDescriptor::IOProcessor();
    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    // Plotfiles.
    const int nPlotFiles = pp.countval("infile");
    if (nPlotFiles <= 0)
      amrex::Abort("No plotfiles specified: set infile=PLT1 [PLT2 ...]");
    Vector<std::string> plotFileNames(nPlotFiles);
    for (int i = 0; i < nPlotFiles; ++i)
      pp.get("infile", plotFileNames[i], i);
    if (verbose)
      std::cout << "Processing " << nPlotFiles << " plotfile(s)." << std::endl;

    // Variables.
    const int nVars = pp.countval("vars");
    if (nVars <= 0)
      amrex::Abort("No variables specified: set vars=\"VAR1 [VAR2 ...]\"");
    Vector<std::string> whichVar(nVars), whichVarOut(nVars);
    for (int v = 0; v < nVars; ++v) {
      pp.get("vars", whichVar[v], v);
      whichVarOut[v] = ProtectSlashes(whichVar[v]);
    }

    // Options.
    int inFinestLevel = -1;
    pp.query("finestLevel", inFinestLevel);

    int goPastMaxI = 1;
    pp.query("goPastMax", goPastMaxI);
    const bool goPastMax = (goPastMaxI != 0);

    int reverseI = 0;
    pp.query("reverse", reverseI);
    const bool reverse = (reverseI != 0);

    std::string format = "ppm";
    pp.query("format", format);
    if (format != "ppm" && format != "png" && format != "pdf")
      amrex::Abort("Unknown format '" + format + "' (use ppm, png or pdf)");
    const std::string ext = "." + format;

    std::string colormap = "jet";
    pp.query("colormap", colormap);
    const bool builtinCmap = (colormap == "jet" || colormap == "grayscale");
    const int cmapIdx =
      builtinCmap ? -1 : pele_colormaps::colormapIndex(colormap);
    if (!builtinCmap && cmapIdx < 0)
      amrex::Abort(
        "Unknown colormap '" + colormap + "'. Available: jet grayscale " +
        pele_colormaps::availableColormaps());

    std::string outputDir = ".";
    pp.query("outputDir", outputDir);
    if (
      ParallelDescriptor::IOProcessor() && !outputDir.empty() &&
      outputDir != ".") {
      if (!UtilCreateDirectory(outputDir, 0755))
        CreateDirectoryFailed(outputDir);
    }
    ParallelDescriptor::Barrier();

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    // Loop over plotfiles.
    for (int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
      const std::string& infile = plotFileNames[iPlot];
      if (verbose)
        std::cout << "\nOpening " << infile << " ..." << std::endl;

      DataServices dataServices(infile, fileType);
      if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
      AmrData& amrData = dataServices.AmrDataRef();

      // Validate variable names (works for stored and derived quantities).
      for (int v = 0; v < nVars; ++v)
        if (!amrData.CanDerive(whichVar[v]))
          amrex::Abort("Variable not in plotfile: " + whichVar[v]);

      // Finest level to use.
      int finestLevel = amrData.FinestLevel();
      if (inFinestLevel > -1 && inFinestLevel < finestLevel)
        finestLevel = inFinestLevel;
      const int nLevels = finestLevel + 1;
      const Vector<Box>& probDomain = amrData.ProbDomain();

      // Color range per variable: data min/max across all levels, unless
      // overridden with useminmax<i>=MIN MAX (i is 1-based).
      Vector<Real> vMin(nVars), vMax(nVars);
      for (int v = 0; v < nVars; ++v) {
        vMin[v] = 1e100;
        vMax[v] = -1e100;
        for (int lev = 0; lev < nLevels; ++lev) {
          Real mn, mx;
          amrData.MinMax(probDomain[lev], whichVar[v], lev, mn, mx);
          vMin[v] = std::min(vMin[v], mn);
          vMax[v] = std::max(vMax[v], mx);
        }
        const std::string argName = "useminmax" + std::to_string(v + 1);
        const int nMinMax = pp.countval(argName.c_str());
        if (nMinMax > 0) {
          if (nMinMax != 2)
            amrex::Abort("Need exactly 2 values for " + argName);
          pp.get(argName.c_str(), vMin[v], 0);
          pp.get(argName.c_str(), vMax[v], 1);
        }
        if (verbose)
          std::cout << "  " << whichVar[v] << " range: " << vMin[v] << " / "
                    << vMax[v] << std::endl;
      }

      // Build the (single) box to render: the whole domain in 2-D, or one
      // slice plane in 3-D.
      Box tempBox(probDomain[finestLevel]);
      int width = 0, height = 0;
#if (AMREX_SPACEDIM == 3)
      std::string sliceTag;
      int islice = -1, sliceDir = -1;
      std::string sliceStr;
      const Box& fineDom = probDomain[finestLevel];
      const char* dirName[3] = {"X", "Y", "Z"};
      const char* dirPfx[3] = {"x", "y", "z"};

      // A slice plane can be given per direction either as a cell index
      // (x/y/zslice) or as a physical coordinate (x/y/zsliceCoord), the latter
      // snapped to the nearest cell of the level in use. At most one direction
      // may be specified.
      for (int d = 0; d < 3; ++d) {
        int idxVal = 0;
        bool haveThis = false;
        const std::string idxOpt = std::string(dirPfx[d]) + "slice";
        const std::string crdOpt = std::string(dirPfx[d]) + "sliceCoord";
        if (pp.countval(idxOpt.c_str()) > 0) {
          pp.get(idxOpt.c_str(), idxVal);
          haveThis = true;
        }
        if (pp.countval(crdOpt.c_str()) > 0) {
          Real coord;
          pp.get(crdOpt.c_str(), coord);
          const Real plo = amrData.ProbLo()[d];
          const Real dx = amrData.DxLevel()[finestLevel][d];
          // Nearest cell center: center of cell j is plo + (j + 0.5) * dx.
          idxVal = (int)std::lround((coord - plo) / dx - 0.5);
          haveThis = true;
        }
        if (haveThis) {
          if (sliceDir != -1)
            amrex::Abort("Specify only one slice plane (one of x/y/zslice or "
                         "x/y/zsliceCoord)");
          sliceDir = d;
          islice = idxVal;
          sliceStr = dirName[d];
        }
      }
      if (sliceDir == -1) { // default: yslice=0
        sliceDir = Amrvis::YDIR;
        islice = 0;
        sliceStr = "Y";
      }
      // Keep the plane inside the domain (guards out-of-range indices and
      // coordinates that round just past the boundary).
      islice = std::max(
        fineDom.smallEnd(sliceDir), std::min(fineDom.bigEnd(sliceDir), islice));
      tempBox.setSmall(sliceDir, islice);
      tempBox.setBig(sliceDir, islice);
      sliceTag = "_" + sliceStr + std::to_string(islice);
      switch (sliceDir) {
      case 0:
        width = probDomain[finestLevel].length(Amrvis::YDIR);
        height = probDomain[finestLevel].length(Amrvis::ZDIR);
        break;
      case 1:
        width = probDomain[finestLevel].length(Amrvis::XDIR);
        height = probDomain[finestLevel].length(Amrvis::ZDIR);
        break;
      default:
        width = probDomain[finestLevel].length(Amrvis::XDIR);
        height = probDomain[finestLevel].length(Amrvis::YDIR);
        break;
      }
      if (verbose)
        std::cout << "  slice " << sliceStr << islice << std::endl;
#else
      width = probDomain[finestLevel].length(Amrvis::XDIR);
      height = probDomain[finestLevel].length(Amrvis::YDIR);
#endif

      // Fill the box on the I/O processor only (keeps image assembly and file
      // writing serial and collective-free).
      BoxArray domainBoxArray(tempBox);
      Vector<int> pmap(1, ParallelDescriptor::IOProcessorNumber());
      DistributionMapping dm(pmap);
      MultiFab mf(domainBoxArray, dm, nVars, 0);

      Vector<int> destFills(nVars);
      for (int v = 0; v < nVars; ++v)
        destFills[v] = v;
      amrData.FillVar(mf, finestLevel, whichVar, destFills);

      // Assemble and write one image per variable (I/O processor only).
      for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const FArrayBox& fab = mf[mfi];
        std::vector<unsigned char> rgb((size_t)width * height * 3);
        for (int v = 0; v < nVars; ++v) {
          const Real* dat = fab.dataPtr(v);
          Real span = vMax[v] - vMin[v];
          if (span <= 0.0)
            span = 1.0;
          const Real cmax =
            (colormap == "jet" && goPastMax && !reverse) ? 1.5 : 1.0;
          for (int k = 0; k < height; ++k) {
            for (int i = 0; i < width; ++i) {
              // Flip vertically so the image is right-side up.
              const size_t px = ((size_t)(height - k - 1) * width + i) * 3;
              Real t = (dat[(size_t)k * width + i] - vMin[v]) / span;
              if (reverse)
                t = 1.0 - t;
              const Real c = std::max(Real(0.0), std::min(cmax, t));
              unsigned char r, g, b;
              if (cmapIdx >= 0)
                pele_colormaps::sampleColormap(cmapIdx, c, r, g, b);
              else if (colormap == "grayscale")
                ColormapGray(c, r, g, b);
              else
                ColormapJet(c, goPastMax, r, g, b);
              rgb[px + 0] = r;
              rgb[px + 1] = g;
              rgb[px + 2] = b;
            }
          }

          std::string fname = PlotfileBasename(infile) + "_" + whichVarOut[v];
#if (AMREX_SPACEDIM == 3)
          fname += sliceTag;
#endif
          fname += ext;
          if (!outputDir.empty() && outputDir != ".")
            fname = outputDir + "/" + fname;

          bool ok = false;
          if (format == "png")
            ok = WritePNG(fname, width, height, rgb);
          else if (format == "pdf")
            ok = WritePDF(fname, width, height, rgb);
          else
            ok = WritePPM(fname, width, height, rgb);
          if (!ok)
            amrex::Abort("Could not write image file: " + fname);
          if (verbose)
            std::cout << "  wrote " << fname << " (" << width << "x" << height
                      << ")" << std::endl;
        }
      }
    }
  }
  amrex::Finalize();
  return 0;
}
