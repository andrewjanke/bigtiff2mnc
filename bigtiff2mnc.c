There are 3 notes (this is a dev code):
1 - TILEAVG is to average each tile in a single voxel. This is for
    debuging only - to have a quick coarse representation of the image.
2 - vymin, vymax, vxmin, vxmax should be the full size of the image.
    I have something in there to make it a subset of the image (crop).
    These images can be very large and too large to load in register/Display
    at full size. For BigBrain at 1-micron, images are 130,000x110,000.
3 - 65535 stuff is to scale the intensities of the BigBrain.

You'll have to clean the code a little bit. :-)

Have fun!

Claude

----------------------------------------------------------------

/*
 * There are 3 notes (this is a dev code):
 * 1 - TILEAVG is to average each tile in a single voxel. This is for
 *     debuging only - to have a quick coarse representation of the image.
 * 2 - vymin, vymax, vxmin, vxmax should be the full size of the image.
 * 
 *     I have something in there to make it a subset of the image (crop).
 *     These images can be very large and too large to load in register/Display
 *     at full size. For BigBrain at 1-micron, images are 130,000x110,000.
 * 
 * 3 - 65535 stuff is to scale the intensities of the BigBrain.
 * 
 * To build
 * 1 - Obtain bigtiff lib package from http://bigtiff.org/#API_CHANGES
 *     Source ZIP for libtiff_4.1 only.
 * 
 * 2 - unzip libtiff-4.1_only.zip
 *     cd tiff-4.1/libtiff
 *     make
 *     cp libtiff.a to minc quarantine Linux-x86_64/lib/
 * 
 * 3 - gcc -DMINC2 -I/home/claude/QUARANTINES/Linux-x86_64/include \
 *     -I/home/claude/QUARANTINES/TIFF/tiff-4.1/libtiff tiff2mnc.c \
 *     -L/home/claude/QUARANTINES/Linux-x86_64/lib -o tiff2mnc -ltiff \
 *     -ljpeg -lminc2 -lnetcdf -lhdf5 -lm -lz
*/


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <tiffio.h>

#include <minc.h>

#define NDIMS 3
#define COMPRESS_LEVEL 4
// #define TILEAVG

extern char *time_stamp(int, char **);

int
main(int argc, char **argv)
{
  TIFF *r_fp;
  int w_fd;
  uint32 x, y, yy, imageWidth, imageLength, tileWidth, tileLength, tilesize;
  uint32 vxmin, vxmax, vymin, vymax;
  uint16 bps;
  double max;
  double min;
  tstrip_t i;
  unsigned short **tiledata;
  unsigned short *mincdata;
  double vrange[2];
  int dimids[NDIMS];
  long count[NDIMS];
  long start[NDIMS];
  int j, ii, jj;
  int var_id;
  char *history_str = time_stamp(argc, argv);
  double xstep, ystep;
  float xres, yres;
  uint16 resunit;

  if (argc < 3) {
    fprintf(stderr, "Not enough arguments\n");
  }

  r_fp = TIFFOpen(argv[1], "r");
  if (r_fp == NULL) {
    fprintf(stderr, "Can't open '%s' for reading\n", argv[1]);
    exit(-1);
  }

  if (TIFFGetField(r_fp, TIFFTAG_IMAGEWIDTH, &imageWidth)) {
    printf("(image width %d) ", imageWidth);
  }

  if (TIFFGetField(r_fp, TIFFTAG_IMAGELENGTH, &imageLength)) {
    printf("(image length %d) ", imageLength);
  }

  if (TIFFGetField(r_fp, TIFFTAG_TILEWIDTH, &tileWidth)) {
    printf("(tile width %d) ", tileWidth);
  }

  if (TIFFGetField(r_fp, TIFFTAG_TILELENGTH, &tileLength)) {
    printf("(tile length %d) ", tileLength);
  }

  if (TIFFGetField(r_fp, TIFFTAG_BITSPERSAMPLE, &bps)) {
    printf("(%d bps) ", bps);
  }

  if (TIFFGetField(r_fp, TIFFTAG_XRESOLUTION, &xres)) {
    // xres = 24000;
    printf("(xres %f) ", xres);
  }

  if (TIFFGetField(r_fp, TIFFTAG_YRESOLUTION, &yres)) {
    // yres = 24000;
    printf("(yres %f) ", yres);
  }

  if (TIFFGetField(r_fp, TIFFTAG_RESOLUTIONUNIT, &resunit)) {
    printf("(unit %d)", resunit);
  }

  printf("\n");

#if 0
  if (resunit == RESUNIT_INCH) {
    printf( "Units are in inches...\n" );
    xstep = 25.4 / xres;        /* convert to mm */
    ystep = 25.4 / yres;        /* convert to mm */
  } else if (resunit == RESUNIT_CENTIMETER) {
    printf( "Units are in cm...\n" );
    xstep = 10.0 / xres;
    ystep = 10.0 / yres;
  }
#else
  xstep = 0.001 * xres;  // units are read in microns, convert to mm
  ystep = 0.001 * yres;
#endif

  // Open the minc file for writing

#ifdef MINC2
  struct mi2opts opts;
  opts.struct_version = MI2_OPTS_V1;
  opts.comp_type = MI2_COMP_ZLIB;
  opts.comp_param = COMPRESS_LEVEL;       // compress image
  opts.chunk_type = MI2_CHUNK_UNKNOWN;    // will use default chunking
  w_fd = micreatex(argv[2], NC_NOCLOBBER|MI2_CREATE_V2, &opts);
#else
  w_fd = micreate(argv[2], NC_NOCLOBBER);
#endif
  if (w_fd < 0) {
    fprintf(stderr, "Can't open '%s' for writing\n", argv[1]);
    exit(-1);
  }

#ifdef TILEAVG
  dimids[0] = ncdimdef(w_fd, MIzspace, 1);
  dimids[1] = ncdimdef(w_fd, MIyspace, (int)(imageLength/tileLength)+1 );
  dimids[2] = ncdimdef(w_fd, MIxspace, (int)(imageWidth/tileWidth)+1 );
#else

#if 0
  vxmin = (int)( ( 0.0 * imageWidth ) / tileWidth ) * tileWidth;
  vxmax = (int)( ( 0.3571 * imageWidth ) / tileWidth ) * tileWidth - 1;
  vymin = (int)( ( 0.0 * imageLength ) / tileLength ) * tileLength;
  vymax = (int)( ( 1.0 * imageLength ) / tileLength ) * tileLength - 1;
#else
  vxmin = 0;
  vxmax = imageWidth-1;
  vymin = 0;
  vymax = imageLength-1;
#endif

  dimids[0] = ncdimdef(w_fd, MIzspace, 1);
  dimids[1] = ncdimdef(w_fd, MIyspace, vymax - vymin + 1 );
  dimids[2] = ncdimdef(w_fd, MIxspace, vxmax - vxmin + 1 );
#endif

  micreate_std_variable(w_fd, MIimage, NC_SHORT, 3, dimids);
  micreate_std_variable(w_fd, MIimagemin, NC_DOUBLE, 0, NULL);
  micreate_std_variable(w_fd, MIimagemax, NC_DOUBLE, 0, NULL);

  micreate_std_variable(w_fd, MIzspace, NC_INT, 0, NULL);
  micreate_std_variable(w_fd, MIyspace, NC_INT, 0, NULL);
  micreate_std_variable(w_fd, MIxspace, NC_INT, 0, NULL);

  miattputstr(w_fd, NC_GLOBAL, MIhistory, history_str);

  var_id = ncvarid(w_fd, MIimage);
  miattputstr(w_fd, var_id, MIsigntype, MI_UNSIGNED);
  miattputstr(w_fd, var_id, MIcomplete, MI_TRUE);

  vrange[0] = 0;
  vrange[1] = USHRT_MAX;

  miset_valid_range(w_fd, var_id, vrange);
  var_id = ncvarid(w_fd, MIxspace);
  miattputstr(w_fd, var_id, MIunits, "mm");
#ifdef TILEAVG
  miattputdbl(w_fd, var_id, MIstep, xstep*tileWidth);
#else
  miattputdbl(w_fd, var_id, MIstep, xstep);
#endif
  miattputdbl(w_fd, var_id, MIstart, -(xstep * imageWidth) / 2.0 );
  printf( "x from %f by %f for %d voxels\n", -(xstep * imageWidth) / 2.0,
          xstep, imageWidth );

  var_id = ncvarid(w_fd, MIyspace);
  miattputstr(w_fd, var_id, MIunits, "mm");
#ifdef TILEAVG
  miattputdbl(w_fd, var_id, MIstep, -ystep*tileLength);
#else
  miattputdbl(w_fd, var_id, MIstep, -ystep);
#endif
  miattputdbl(w_fd, var_id, MIstart, (ystep * imageLength) / 2.0 );
  printf( "y from %f by %f for %d voxels\n", (ystep * imageLength) / 2.0,
          -ystep, imageLength );

  var_id = ncvarid(w_fd, MIzspace);
  miattputstr(w_fd, var_id, MIunits, "mm");
  miattputdbl(w_fd, var_id, MIstep, 0.02);
  miattputdbl(w_fd, var_id, MIstart, 0.0);
  printf( "z from %f by %f for %d voxels\n", 0.0, 0.02, 1 );

  ncendef(w_fd);

  max = 0;
  min = USHRT_MAX;

  var_id = ncvarid(w_fd, MIimage);

  count[0] = 1;  // number of slices per tif
  count[1] = 1;  // number of rows per strip
#ifdef TILEAVG
  count[2] = (int)(imageWidth/tileWidth);  // number of tiles per row
#else
  count[2] = vxmax - vxmin + 1;  // number of columns per row
#endif
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  tilesize = TIFFTileSize(r_fp);
  printf("tilesize %d\n", tilesize);

  short numtiles = ( imageWidth + tileWidth - 1 ) / tileWidth;
  printf( "numtiles = %d\n", numtiles );
  tiledata = (unsigned short**)malloc( numtiles * sizeof( unsigned short * ) );
  for( j = 0; j < numtiles; j++ ) {
    tiledata[j] = (unsigned short *)malloc( tilesize );
    if( tiledata[j] == NULL ) {
      printf( "Cannot allocate memory for tile %d\n", j );
      exit(1);
    }
  }

  mincdata = (unsigned short * )malloc( imageWidth * sizeof( unsigned short ));

  for( y = 0; y < imageLength; y += tileLength ) {
    printf( "Processing tiles at y = %d...\n", y );
    for( x = 0; x < imageWidth; x += tileWidth ) {
      TIFFReadTile( r_fp, tiledata[x/tileWidth], x, y, 0, 0 );
    }
#ifdef TILEAVG
    // average each tile into a voxel
    for( x = 0; x < imageWidth; x += tileWidth ) {
      double avg = 0.0;
      j = x/tileWidth;
      for( ii = 0; ii < tileLength; ii++ ) {
        for( jj = 0; jj < tileWidth; jj++ ) {
          unsigned short val = tiledata[j][ii*tileWidth+jj];
          if (bps == 8) {
            /* scale up from byte to short */
            unsigned char b = (unsigned char)val;
            val = (unsigned short)(((double)b / 255.0)* 65535.0);
          }
          avg += val;
        }
      }
      avg = ( avg / tileLength ) / tileWidth;
      if (avg > max) {
        max = avg;
      }
      if (avg < min) {
        min = avg;
      }
      mincdata[j] = 65535.0 - avg;
    }
    mivarput( w_fd, var_id, start, count, NC_SHORT, MI_UNSIGNED, mincdata );
    start[1] += count[1];
#else
    // read line-by-line from the tiles
    for( yy = 0; yy < tileLength && y + yy < imageLength; yy++ ) {
      if( y+yy >= vymin && y+yy <= vymax ) {
        for( x = 0; x < imageWidth; x++ ) {
          unsigned short val = tiledata[x/tileWidth][yy*tileWidth+(x%tileWidth)];
          if (bps == 8) {
            /* scale up from byte to short */
            unsigned char b = (unsigned char)val;
            val = (unsigned short)(((double)b / 255.0)* 65535.0);
          }
          if (val > max) {
            max = val;
          }
          if (val < min) {
            min = val;
          }
          mincdata[x] = 65535.0 - val;
        }
        mivarput( w_fd, var_id, start, count, NC_SHORT, MI_UNSIGNED, &mincdata[vxmin] );
        start[1] += count[1];
      }
    }
#endif
  }
  TIFFClose(r_fp);

  for( j = 0; j < numtiles; j++ ) {
    free( tiledata[j] );
  }
  free( tiledata );
  free( mincdata );

  start[0] = 0;
  mivarput1(w_fd, ncvarid(w_fd, MIimagemin), start, NC_DOUBLE, MI_SIGNED,
            &min);
  start[0] = 0;
  mivarput1(w_fd, ncvarid(w_fd, MIimagemax), start, NC_DOUBLE, MI_SIGNED,
            &max);

  miclose(w_fd);
}
