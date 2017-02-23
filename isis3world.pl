#!/usr/bin/perl -s
###############################################################################
#
#_TITLE  isis3world.pl - Create a GIS detached header and/or GIS world file
#                           for ISIS 3 image. ISIS 3 cubes should be exported to
#                           a BSQ format or other format like TIFF or Jpeg.
#                           The ISIS 3 cube will most
#                           likely also need a new extension. This program
#                           does not change the extension for you.
#
# Note this version does not need ISIS 3 installed to run
#
#  Pieces of code used are originally from Dr. Herrick
#
#_ARGS  
#
#_Parm  Output header bit types ( 8 or 16)
#
#_Parm  Output header/world formats - gxp (mhdr), raw( hdr), tfw, gfw, jgw, pgw, cube(hdr), aux
#
#_Parm  Input filename
#
#_USER  Command line entry options [optional parameters]:
#
#   isis3world.pl [-bit=8|16] [-r|-g|-t|-c|-j|-p|-P|-gxp] input.cub
#
# Requirements:
#
#   - images should be level 2 (map projected)
#
# Examples: 
#       1) To generate a tif worldfile 
#
#              isis3world.pl -t input.cub
#
#       2) To generate an 8-bit ISIS cube file using all defaults:
#              isis3world.pl input.bsq
#
#_DESC Convert an ISIS 3 cube to a GIS header and world files.
#
#      Options:
#       1) Input file from 8/16 bit (-bit=8|16)
#             Default: 8 bit
#
#       2) Output Options:
#                (-r) Generate a raw output header and world files
#
#                (-gxp) Generate a raw output header for Socet GXP
#
#                (-e) Generate a ERDAS raw output header and world files
#
#                (-g) Generate gif wold file
#
#                (-t) Generate tif world file
#
#                (-j) Generate jpeg world file
#
#                (-J) Generate jpeg2000 world file
#
#                (-P) Generate png world file
#
#                (-p) Generate a PCI Aux header (w/ georef) file
#
#                (-prj) Generate a ESRI Well Known Text Projection file (useful for ArcMap and GDAL)
#
#                (-o=input.raw,8|16) an override swtich ISIS file to point at raw file instead of an ISIS file.
#
#             Default:
#                (-c) Generate an ISIS output cube header (w/georef) file
#                
#       3) Output GIS header filenames:  input.hdr, input_8b.hdr, input_16b.hdr, input.raw, input.aux
#                      world GIS files:  input.rrw, input.tfw, input.gfw, input.jgw, input.pgw
#
#_FILE Input files utilized 
#
#
#_LIMS  Design limitations
#        - General use
#    
#_CALL  List of calls:
#
#       List of external Perl modules
#          Math::Trig
#
#       List of ISIS programs that are called:
#
#       Unix utility: 
#
#_HIST
#       Jan 07 2005 - Trent Hare - Program rewritten from isis2world.pl
#       Jun 19 2006 - T.H. - fix PCI Aux bug and add png support
#       Feb 14 2007 - T.H. - fix skipbyte error and another PCI Aux bug
#       Sep 16 2007 - T.H. - Added "prj" option to create a ESRI projection file
#                          - Added "-J" for Jpeg2000 *.j2w
#       Sep 17 2007 - T.H. - Added "-o=input.raw,8|16" to override output with raw
#       Oct 25 2007 - T.H. - Added ERDAS raw TILED output
#       Nov 29 2007 - T.H. - Added support for uppercase projections from ISIS3 label
#       Nov 30 2007 - T.H. - Added back in support for PCI Aux to simulate a 2 band
#                            image as a 3 band for anaglyph conversion
#       July 2 2009 - T.H. - Change Standard_Parallel_1 to latitude_of_origin for Polar Stereographic
#       Mar 25 2015 - T.H. - Added GXP mhdr format
#       May 18 2016 - T.H. - added support for multiple extensions
#       Feb 23 2017 - T.H. - Updated Simple Cylindrical to force a sphere since that is what ISIS uses in
#                                                the projection (for either ographic or ocentric).

#FORMAT TILED
#TILE WIDTH 64
#TILE HEIGHT 64

#
#_END
######################################################################
#includes
use Math::Trig;

######################################################################
# For help - user can enter this perl script and return
######################################################################
   if ($#ARGV < 0)
      {
      print " \n\n          *** HELP ***\n\n";
      print "isis3world.pl -  Create GIS header and world files from an ISIS 3 cube\n\n";
      print "Command line: \n";
      print "  isis3world.pl [-bit=8|16] [-r|-e|-g|-t|-c|-ji|-J|-p|-P|-gxp] [-o=input.raw,8|16] input.cub\n";
      print "    -r = output raw header w/ georefencing (8, 16 bit)\n";
      print "    -gxp = output Socet GXP header no georefencing and GIS worldfile (not used)\n";
      print "    -e = output ERDAS raw header and world file (8, 16, 32 bit)\n";
      print "    -g = output gif world file\n";
      print "    -t = output tif world file\n";
      print "    -J = output jpeg2000 world file\n";
      print "    -j = output jpeg world file\n";
      print "    -P = output png world file\n";
      print "    -p = output PCI Aux header w/ georefencing (8, 16, 32 bit)\n";
      print "    -p=rgg output PCI Aux header repeating last band.  Use to make 2 band anaglyph to look like 3 bands\n";
      print "    -c = output cub header w/ georef (default) (8, 16 bit)\n";
      print "\n    -prj = create ESRI Well Known Text projection file *.prj\n";
      print "\n    -o = Override input ISIS3 file with raw file and change to 0 skipbytes\n";

      print "\nExamples:\n";
      print "   Create header for SocetGXP: isis3world.pl -gxp input.cub\n";
      print "   Create header for 16 bit cube: isis3world.pl -bit=16 -r input.cub\n";
      print "   Create files for 32 bit ERDAS: isis3world.pl -e input.cub\n";
      print "   Create header for 32 bit PCI Aux: isis3world.pl -p input.cub\n";
      print "   Create header for 8 bit PCI Aux with input file override: isis3world.pl -p -o=input.raw,8 input.cub\n";
      print "   Create files for tif: isis3world.pl -t input.cub\n";
      print "   Apply all defaults:     isis3world.pl input.cub\n\n";
      exit 1;
      }

   $input = $ARGV[0];
   chomp $input;

######################################################################
#  Check input file name for .cub extention
######################################################################

   #@fname = split('\.',$input);
   #$root = $fname[0];
   (my $root = $input) =~ s/\.[^.]+$//;
   #print "root";

   #  Create the header file from the bil file
   #print "$input";
   open INIMAGE, $input;

######################################################################
#  Check for filetype if none then default to cube
######################################################################

   if (!($t) && !($g) && !($j) && !($r) && !($c) && !($e) && !($p) && !($P) && !($J) && !($gxp)) {
     $c = 1; 
   }


######################################################################
#  Check for system     - Not needed for this
######################################################################
#
#   $thesys = `uname`;
#   chomp $thesys;          #SunOS or Linux fo ISIS compatible systems
#
##################################################################
#  If the user set the -bit flag, double check the entered values
#     8 and 16 are the only valid values.
#  If -bit was not flagged, then assume an 8bit file as input
##################################################################

   if ($bit) 
    {
     if ($bit !=8 && $bit != 16)
       {
        print "[isis3world.pl-ERROR] Bit switch must be set to 8 or 16\n";
        print " Example:  isis3world.pl -r -bit=16 input.cub\n";
        exit;
       }
    }

##################################################################
# Set defaults if not Image Map is found
##################################################################
$xdim = 1;
$ydim = -1;
$sp_offset = 0.5;
$lp_offset = 0.5;
$tileSamples = 0;
$tileLines = 0;

##################################################################
#Set ISIS NODATA
##################################################################
$NULL1=0;
$NULL2=-32768;
$NULL4=-0.3402822655089E+39;
$tslat = 0; #Set trueScaleLat to 0 as the default, some older ISIS files may not have it

##################################################################
# Loop through image header to extract info
##################################################################
$flag = 0;
while (<INIMAGE>) {

   ###############################################################
   # Extract the SKIPBYTES from Record_bytes and ^Qube
   ###############################################################
    if (/ StartByte /) {
       @sb = split(/ = /,$_);
       $skipbytes = @sb[1];
       chomp $skipbytes;
       $skipbytes = $skipbytes -1;
    }

   ###############################################################
   # Extract the FORMAT style
   ###############################################################
    if (/ Format /) {
        @ft = split(/ = /,$_) ;
        $format = $ft[1];
        chomp $format;
        ###############################################################
        # Warning if the the cube is Tiled. This will not work for external software. You have to convert to BSQ.
        ###############################################################
        if (($r) || ($c)  || ($p)) {
            if (($format eq "Tile") && (!($o))) {
              print "[isis3world.pl-WARNING] You are writing a header for a Tiled ISIS image. Convert ISIS image to BSQ first.\n";
            }
        }
    }

    if (/ TileSamples /) {
       @t1 = split(/ = /,$_) ;
       $tileSamples = $t1[1];
       chomp $tileSamples;
    }

    if (/ TileLines /) {
       @t1 = split(/ = /,$_) ;
       $tileLines = $t1[1];
       chomp $tileLines;
    }


    ###############################################################
    # Extract the sample and line dimension of input file
    ###############################################################
    if (/ Samples /) {
       @t1 = split(/ = /,$_) ;
       $samples = $t1[1];
       chomp $samples;
    }
    if (/ Lines /) {
       @t1 = split(/ = /,$_) ;
       $lines = $t1[1];
       chomp $lines;
    }
        if (/ Bands /) {
       @t1 = split(/ = /,$_) ;
       $bands = $t1[1];
       chomp $bands;
    }
    
     ###############################################################
     # Extract the CORE_ITEM_BYTES
     ###############################################################
     if (/  Type  /) {
        if ($flag == 0) { #This is a hacky fix to not catch "Type" again.
          @itp = split(/ = /,$_) ;
          $it = $itp[1];
          chomp $it;
          if ($it eq "UnsignedByte")  {
            $itype = 8;
            #Options U1,U2,U4,U8,U16,S16,U32,F32,F64
            $ERDASitype = "U8";
            #Options 8u,16U,16S,32R
            $PCIitype = "8U";
            $nodata = $NULL1;
            #Remap for PCI Aux and GXP
            $it = 1;
          }
          elsif ($it eq "SignedWord") {
            $itype = 16;
            $ERDASitype = "S16";
            $PCIitype = "16S";
            $nodata = $NULL2;
            #Remap for PCI Aux and GXP
            $it = 2;
          }
          elsif ($it eq "Real") {
            $itype = 32;
            $ERDASitype = "F32";
            $PCIitype = "32R";
            $nodata = $NULL4;
            #Remap for PCI Aux and GXP
            $it = 4;
          }
          else { 
            $itype = 0; 
            $nodata = 0;
            #Remap for PCI Aux
            $it = 0;
          }
       }
       $flag=1;
     }

     ##############################################################
     # Extract the machine type
     ##############################################################
     if (/ ByteOrder /) {
        @tp = split(/ = /,$_) ;
        $type = $tp[1];
        chomp $type;
        if ($type eq "Lsb") {
            $machine = "I";
            $ERDASmachine = "LSB";
            $PCImachine = "Swapped";
            $GXPmachine = "0";
        }
        else {
            $machine = "M";
            $ERDASmachine = "MSB";
            $PCImachine = "Unswapped";
            $GXPmachine = "1";
        }
     }

     ##############################################################
     # Extract the MAP_SCALE
     ###############################################################
     if (/ PixelResolution /) {
         @xydim= split(/ = /,$_) ;
         $xdim = $xydim[1] * 1;
         $ydim = $xydim[1] * -1;
     }
   
     ##############################################################
     # Extract the SAMPLE_PROJECTION_OFFSET / LINE_PROJECTION_OFFSET
     ##############################################################

     if (/ UpperLeftCornerX /) {
        @ulxmap = split(/ = /,$_) ;
        $ulcornerxmap = @ulxmap[1];
	$ulcornerxmap = $ulcornerxmap * 1; #poor man conversion to number
      }
      if (/ UpperLeftCornerY /) {
         @ulymap = split(/ = /,$_) ;
         $ulcornerymap = @ulymap[1];
 	 $ulcornerymap = $ulcornerymap * 1; #poor man conversion to number
       }

     ##############################################################
     # Extract the Projection Keywords
     ##############################################################
     if (/EquatorialRadius/) {
        @a_axis = split(/ = /,$_) ;
        $a_axis = @a_axis[1] * 1;
        #$a_axis = $a_axis.replace("<meters>") * 1;
     }
     if (/PolarRadius/) {
        @c_axis = split(/ = /,$_) ;
        $c_axis = @c_axis[1] * 1;
        #$c_axis = $c_axis.replace("<meters>") * 1;
     }
     if (/CenterLongitude/) {
        @clon = split(/ = /,$_) ;
        $clon = @clon[1] * 1;
     }
     if (/CenterLatitude /) {
        @clat = split(/ = /,$_) ;
        $clat = @clat[1] * 1;
     }
     if (/TrueScaleLatitude/) {
        @tslat = split(/ = /,$_) ;
        $tslat = @tslat[1] * 1;
     }
     if (/TargetName/) {
        @target = split(/ = /,$_) ;
        $target = @target[1];
        chomp $target;
        $target =~ s/\"//;
        $target =~ s/\"//;
     }
     if (/LongitudeDirection/) {
        @londir = split(/ = /,$_) ;
        $londir = @londir[1];
        chomp $londir;
        $londir =~ s/\"//;
        $londir =~ s/\"//;
     }
     if (/LongitudeDomain/) {
        @lonsys = split(/ = /,$_) ;
        $lonsys = @lonsys[1] * 1;
     }
     if (/LatitudeType/) {
        @latType = split(/ = /,$_) ;
        $latType = @latType[1];
        chomp $latType;
        $latType =~ s/\"//;
        $latType =~ s/\"//;
     }
     if (/FirstStandardParallel/) {
        @firstParallel = split(/ = /,$_) ;
        $firstParallel = @firstParallel[1] * 1;
     }
     if (/SecondStandardParallel/) {
        @secondParallel = split(/ = /,$_) ;
        $secondParallel = @secondParallel[1] * 1;
     }
     if (/ScaleFactor/) {
        @scaleFactor = split(/ = /,$_) ;
        $scaleFactor = @scaleFactor[1] * 1;
     }
     if (/ProjectionName/) {
        @projection = split(/ = /,$_) ;
        $projection = @projection[1];
        chomp $projection;
        $projection =~ s/\"//;
        $projection =~ s/\"//;
        #print $projection;
     }
     
     ##############################################################
     # Find the End
     ##############################################################
     if (/Object = Label/) {
        last;
     }

# End of image header loop
}

################################################################
# Do map info prep
################################################################
#if -bit is sent, it will override what the input cub actually is.
if ($bit) {
  $itype = $bit;
}

#Calculate from upper left to center of pixel in meters
$ulxmap = $ulcornerxmap + ($xdim / 2);
$ulymap = $ulcornerymap + ($ydim / 2);

#Check this to make sure ISIS3 is doing what I think it is doing
if (lc($projection) eq "equirectangular") {
     if ($latType eq "Planetocentric") {
         $tslat = $clat;
      }
}

#Calculate localRadius using ISIS3 simple elliptical method 
#  not more standard Radius of Curvature method
if ($tslat != 0) { 
    #calculate new local radius
    $radlat = $tslat * pi / 180;

    $localradius = $a_axis * $c_axis / (sqrt((($c_axis*cos($radlat))**2) + (($a_axis*sin($radlat))**2))); 
    print "Calculated new local radius: $localradius, old radius: $a_axis\n";
    $a_axis = $localradius; 
    $c_axis = $localradius; 
}         

####################################################################
# Swtich input to raw filename and change skipbytes to 0
####################################################################
#if -o is sent, it will override the input file actually is.
#  This helps to stretch an image in ISIS but still read ISIS label for projection
if ($o) {
  @rawInput = split(',',$o);
  $input = $rawInput[0];
  chomp $input;
  $itype = $rawInput[1]  * 1;
  if  (!(-e $input)) {
    print "[isis3world.pl - ERROR] The override input raw file does not exist: $o\nExiting\n";
    exit 1;
  }
  print "[isis3world.pl-WARNING] An override to the input file was sent. file=$input; type=$itype\n";
  $skipbytes = 0;

  if ($itype == 8)  {
    $ERDASitype = "U8";
    $PCIitype = "8U";
    $nodata = $NULL1;
    #Remap for PCI Aux and GXP
    $it = 1;
  }
  elsif ($itype == 16) {
    $ERDASitype = "S16";
    $PCIitype = "16S";
    $nodata = $NULL2;
    #Remap for PCI Aux and GXP
    $it = 2;
  }
  elsif ($itype == 32) {
    $ERDASitype = "F32";
    $PCIitype = "32R";
    $nodata = $NULL4;
    #Remap for PCI Aux and GXP
    $it = 4;
  }
  else {
    print "[isis3world.pl - ERROR] The override bit type must be 8, 16 or 32.\nExiting\n";
    exit 1;
  }
}


####################################################################
# Generate raw or cube file ESRI header (w/ georef) & No world file
####################################################################
  if (($r) || ($c))
  {
    ###############################################
    # define world and header files 
    ###############################################
    #define WORLD file name
    if ($r) {
      #if ($bands > 1) {
      #  $outworld = $root.".bqw";
      #}
      #else {
      #  $outworld = $root.".blw";
      #}
      $outhdr = $root.".hdr";

    } else {
      if ($bit) { #if cube is requested and bit sent, then add obit to name. 
        if ($bit == 16) {
            $obit = "_16b";
        }
        else {
            $obit = "_8b";
        }
        if ($bands > 1) {
          #$outworld = $root.$obit.".bqw";
          $outhdr = $root.$obit.".hdr";
        }
        else {
          #$outworld = $root.$obit.".blw";
          $outhdr = $root.$obit.".hdr";
        }
      } else { #if cube is requested and no bit sent, then add use basename. 
        if ($bands > 1) {
          #$outworld = $root.".bqw";
          $outhdr = $root.".hdr";
        }
        else {
          #$outworld = $root.".blw";
          $outhdr = $root.".hdr";
        }
      }
      if ($bands > 1) {
        print "\n **NOTE: Rename your ISIS cub file to have an extension '.bsq'";
      } else {
        print "\n **NOTE: Rename your ISIS cub file to have an extension '.bil'";
      }
    }

    ###############################################
    # write out header 
    ###############################################
    unlink ("$outhdr");

    open OUTHDR, "> $outhdr ";
    print OUTHDR "nrows $lines\n";
    print OUTHDR "ncols $samples\n";
    print OUTHDR "nbands $bands\n";
    print OUTHDR "nbits $itype\n"; 
    print OUTHDR "byteorder $machine\n";

    if ($bands > 1) {
      print OUTHDR "layout BSQ\n";
    } else {
      print OUTHDR "layout BIL\n";
    }

    if ($r) {
      print OUTHDR "skipbytes 0\n";
    } else {
      print OUTHDR "skipbytes $skipbytes\n";
    }
    print OUTHDR "ulxmap $ulxmap \n";
    print OUTHDR "ulymap $ulymap \n";
    print OUTHDR "xdim $xdim\n";
    $ydim2 = $ydim * -1;
    print OUTHDR "ydim $ydim2\n";
    print OUTHDR "nodata $nodata \n";

    close OUTHDR;


    if (-e $outhdr)
      {
        print "\n Header file generated (no world file needed): $outhdr\n";
      }
    else
     {
       print "\n Header file not generated...something's wrong\n";
       exit;
     }
    }
  
####################################################################
# Generate raw ERDAS header and world
####################################################################
  if (($e))
  {
    ###############################################
    # define world and header files 
    ###############################################
    #define WORLD file name
    $outworld = $root.".rww";
    $outhdr = $root.".raw";



    ###############################################
    # write out header 
    ###############################################
    unlink ("$outhdr");

    open OUTHDR, "> $outhdr ";
    print OUTHDR "IMAGINE_RAW_FILE\n";
    print OUTHDR "WIDTH $samples\n";
    print OUTHDR "HEIGHT $lines\n";
    print OUTHDR "NUM_LAYERS $bands\n";
    if ($format eq "Tile") {
       print OUTHDR "FORMAT TILED\n";
       print OUTHDR "TILE WIDTH $tileSamples\n";
       print OUTHDR "TILE HEIGHT $tileLines\n";
    } else {
       print OUTHDR "FORMAT BSQ\n";
    }
    print OUTHDR "DATA_TYPE $ERDASitype\n"; 
    print OUTHDR "BYTE_ORDER $ERDASmachine\n";
    print OUTHDR "PIXEL_FILES $input\n";
    print OUTHDR "DATA_OFFSET $skipbytes\n";
    close OUTHDR;


    if (-e $outhdr)
      {
        print "\n Header file generated: $outhdr\n";
      }
    else
     {
       print "\n Header file not generated...something's wrong\n";
       exit;
     }
   }

####################################################################
# Generate raw SOCET header (mhdr)
####################################################################
  if (($gxp))
  {
    ###############################################
    # define header file 
    ###############################################
    $outhdr = $root.".mhdr";
    #define WORLD file name 
    #skipping for now. Information defined in keywords.lis (from ISIS3).
    #$outworld = $root.".wld";


    ###############################################
    # write out header 
    ###############################################
    unlink ("$outhdr");

    open OUTHDR, "> $outhdr ";
    print OUTHDR "SOCET_GXP\n";
    print OUTHDR "samples = $samples\n";
    print OUTHDR "lines = $lines\n";
    print OUTHDR "header offset = $skipbytes\n";
    print OUTHDR "file type = ISIS3\n";
    print OUTHDR "data type = $it\n"; 
    print OUTHDR "interleave = BSQ\n";
    print OUTHDR "byte order = $GXPmachine\n";
    if ($format eq "Tile") {
       print OUTHDR "tile lines = $tileLines\n";
       print OUTHDR "tile samples = $tileSamples\n";
    } else {
       print OUTHDR "tile lines = 1\n";
       print OUTHDR "tile samples = 1\n";
    }
    print OUTHDR "default bands = {}\n";
    print OUTHDR "band names = {Unknown}\n";
    print OUTHDR "bands = $bands\n\n";
    close OUTHDR;


    if (-e $outhdr)
      {
        print "\n Header file generated: $outhdr\n";
      }
    else
     {
       print "\n Header file not generated...something's wrong\n";
       exit;
     }
   }

####################################################################
# Generate raw PCI Aux header 
####################################################################
  if (($p))
  {
    ###############################################
    # define header file (w/ georef) 
    ###############################################
    $outhdr = $root.".aux";

    ###############################################
    # write out header 
    ###############################################
    unlink ("$outhdr");

    #Calculate lower right corner from upper right corner of pixel in meters
    $lrcornerxmap = $ulcornerxmap + ($xdim * ($samples));
    $lrcornerymap = $ulcornerymap - ($xdim * ($lines));

    open OUTHDR, "> $outhdr ";
    print OUTHDR "AuxilaryTarget: $input\n";

    #Create a false 3 band image from 2 bands by repeating last band
    if ($p eq "rgg") { 
      $bands_temp = $bands + 1;
      print OUTHDR "RawDefinition: $samples $lines $bands_temp\n";
    }
    else {
      print OUTHDR "RawDefinition: $samples $lines $bands\n";
    }
    for ($i=1; $i<=$bands; $i++) {
      $skip = ($samples * $lines * $it * ($i - 1)) + $skipbytes;
      $lskip = $samples * $it;
      print OUTHDR "ChanDefinition-$i: $PCIitype $skip $it $lskip $PCImachine\n";
    }
    if ($p eq "rgg") { 
      print OUTHDR "ChanDefinition-$i: $PCIitype $skip $it $lskip $PCImachine\n";
    }
    print OUTHDR "UpLeftX: $ulcornerxmap\n";
    print OUTHDR "UpLeftY: $ulcornerymap\n";
    print OUTHDR "LoRightX: $lrcornerxmap\n"; 
    print OUTHDR "LoRightY: $lrcornerymap\n";
    for ($i=1; $i<=$bands; $i++) {
      print OUTHDR "METADATA_IMG_$i\_NO_DATA_VALUE: $nodata \n";
    }
    if ($p eq "rgg") { 
      print OUTHDR "METADATA_IMG_$i\_NO_DATA_VALUE: $nodata \n";
    }
    close OUTHDR;

    if (-e $outhdr)
      {
        print "\n Header file generated: $outhdr\n";
      }
    else
     {
       print "\n Header file not generated...something's wrong\n";
       exit;
     }
    }
  

####################################################################
# Generate for tiff 
####################################################################
  if ($t)
    {
      $outworld = $root.".tfw";
    }
      
####################################################################
# Generate for jpeg 
####################################################################
  if ($j)
    {
      $outworld = $root.".jgw";
    }
####################################################################
# Generate for jpeg2000 
####################################################################
  if ($J)
    {
      $outworld = $root.".j2w";
    }
####################################################################
# Generate for gif 
####################################################################
  if ($g)
    {
      $outworld = $root.".gfw";
    }
####################################################################
# Generate for png 
####################################################################
  if ($P)
    {
      $outworld = $root.".pgw";
    }

####################################################################
# Write out the World file 
####################################################################
   if ($outworld) {
     unlink ("$outworld");
     open OUTWORLD, "> $outworld ";
     print OUTWORLD "$xdim\n";
     print OUTWORLD "0.0\n";
     print OUTWORLD "0.0\n";
     print OUTWORLD "$ydim\n";
     print OUTWORLD "$ulxmap \n";
     print OUTWORLD "$ulymap \n";
     close OUTWORLD;
     if (-e $outworld)
        {
          print " World file generated: $outworld\n\n";
        }
      else
       {
         print " World file not generated...something's wrong\n\n";
         exit;
        }
    }
##################################################################
#  Find out if the user set the -prj projection flag
##################################################################
     # Projection types
     #############################################################
    ## class   Isis::Equirectangular 
    ## class   Isis::LambertConformal 
    ## class   Isis::Mercator 
    ## class   Isis::ObliqueCylindrical 
    ## class   Isis::Orthographic 
    ## class   Isis::PolarStereographic 
    ## class   Isis::SimpleCylindrical 
    ## class   Isis::Sinusoidal 
    ## class   Isis::TransverseMercator 
   if ($prj) 
   {
        ###############################################
        # define projection file 
        ###############################################
        $outprj = $root.".prj";
        unlink ("$outprj");
        open OUTPROJ, "> $outprj ";
        
        $flattening = ($a_axis - $c_axis)/$a_axis;
        if ($flattening != 0) {
            $flattening = 1/$flattening;
        }
        
        if ($londir eq "PositiveWest") { #swap from west to east
          if ($clon > 0) {
            $clon = $clon - 180;
          } else {
            $clon = $clon * -1;
          }  
        }

        #print $target;
        if (lc($projection) eq "sinusoidal") {
            $projection = "Sinusoidal";
            #ISIS uses a spherical equation for this projection so force a sphere
            print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."_Sphere\",DATUM[\"D_".$target."_Sphere\",SPHEROID[\"".$target."_Sphere\",$a_axis,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],UNIT[\"Meter\",1.0]]";
        }   
        elsif (lc($projection) eq "polarstereographic") {
            $projection = "Polar_Stereographic";
            if (($latType eq "None") || ($latType eq "Planetographic")) {
              print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."\",DATUM[\"D_".$target."\",SPHEROID[\"".$target."\",$a_axis,$flattening]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"latitude_of_origin\",$clat],UNIT[\"Meter\",1.0]]";
            } else {
              #for ocentric latitude system use polar radius as a sphere
              print OUTPROJ "PROJCS[\"".$target."_".$projection."_Sphere\",GEOGCS[\"GCS_".$target."_Sphere_Polar\",DATUM[\"D_".$target."_Sphere_Polar\",SPHEROID[\"".$target."_Sphere_Polar\",$c_axis,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"latitude_of_origin\",$clat],UNIT[\"Meter\",1.0]]";
            }
        }   
        elsif (lc($projection) eq "simplecylindrical") {
            $projection = "Equidistant_Cylindrical";
            #ISIS uses a spherical equation for this projection so force a sphere
            print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."\",DATUM[\"D_".$target."\",SPHEROID[\"".$target."\",$a_axis,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Standard_Parallel_1\",0.0],UNIT[\"Meter\",1.0]]";
        }   
        elsif (lc($projection) eq "equirectangular") {
            $projection = "Equidistant_Cylindrical";
            if (($latType eq "None") || ($latType eq "Planetographic")) {
              print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."\",DATUM[\"D_".$target."\",SPHEROID[\"".$target."\",$a_axis,$flattening]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Standard_Parallel_1\",$clat],UNIT[\"Meter\",1.0]]";
            } else {
              #for ocentric latitude system use semi-major radius as a sphere. Note local radius was calculated above for this projection.
              print OUTPROJ "PROJCS[\"".$target."_".$projection."_Sphere\",GEOGCS[\"GCS_".$target."_Sphere\",DATUM[\"D_".$target."_Sphere\",SPHEROID[\"".$target."_Sphere\",$a_axis,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Standard_Parallel_1\",$clat],UNIT[\"Meter\",1.0]]";
            }
        }   
        elsif (lc($projection) eq "orthographic") {
            $projection = "Orthographic";
            #ISIS uses a spherical equation for this projection so force a sphere
            print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."_Sphere\",DATUM[\"D_".$target."_Sphere\",SPHEROID[\"".$target."_Sphere\",$a_axis,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],UNIT[\"Meter\",1.0]]";
        }   
        elsif (lc($projection) eq "stereographic") {
            $projection = "Stereographic";
            #ISIS uses a spherical equation for this projection so force a sphere
            print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."\",DATUM[\"D_".$target."\",SPHEROID[\"".$target."\",$a_axis,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Latitude_of_Origin\",0.0],PARAMETER[\"Scale_Factor\",1.0],UNIT[\"Meter\",1.0]]";
        }   
        elsif (lc($projection) eq "transversemercator") {
            $projection = "Transverse_Mercator";
            if (($latType eq "None") || ($latType eq "Planetographic")) {
              print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."_Sphere\",DATUM[\"D_".$target."_Sphere\",SPHEROID[\"".$target."_Sphere\",$a_axis,$flattening]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Scale_Factor\",$scaleFactor],PARAMETER[\"Latitude_Of_Origin\",$clat],UNIT[\"Meter\",1.0]]";
            } else {
              #for ocentric latitude system use semi-major radius as a sphere
              print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."_Sphere\",DATUM[\"D_".$target."_Sphere\",SPHEROID[\"".$target."_Sphere\",$a_axis,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Scale_Factor\",$scaleFactor],PARAMETER[\"Latitude_Of_Origin\",$clat],UNIT[\"Meter\",1.0]]";
            }
        }   
        elsif (lc($projection) eq "mercator") {
            $projection = "Mercator";
            if (($latType eq "None") || ($latType eq "Planetographic")) {
              print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."\",DATUM[\"D_".$target."\",SPHEROID[\"".$target."\",$a_axis,$flattening]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Standard_Parallel_1\",$clat],UNIT[\"Meter\",1.0]]";
            } else {
              #for ocentric latitude system use semi-major radius as a sphere
              print OUTPROJ "PROJCS[\"".$target."_".$projection."_Sphere\",GEOGCS[\"GCS_".$target."_Sphere\",DATUM[\"D_".$target."_Sphere\",SPHEROID[\"".$target."_Sphere\",$a_axis,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Standard_Parallel_1\",$clat],UNIT[\"Meter\",1.0]]";
            }
        }   
        elsif (lc($projection) eq "lambertconformal") {
            $projection = "Lambert_Conformal_Conic";
            if (($latType eq "None") || ($latType eq "Planetographic")) {
              print OUTPROJ "PROJCS[\"".$target."_".$projection."\",GEOGCS[\"GCS_".$target."\",DATUM[\"D_".$target."\",SPHEROID[\"".$target."\",$a_axis,$flattening]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Standard_Parallel_1\",$firstParallel],PARAMETER[\"Standard_Parallel_2\",$secondParallel],PARAMETER[\"Latitude_Of_Origin\",$clat],UNIT[\"Meter\",1.0]]";
            } else {
              #for ocentric latitude system use semi-major radius as a sphere
              print OUTPROJ "PROJCS[\"".$target."_".$projection."_Sphere\",GEOGCS[\"GCS_".$target."_Sphere\",DATUM[\"D_".$target."_Sphere\",SPHEROID[\"".$target."_Sphere\",$a_axis,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"".$projection."\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",$clon],PARAMETER[\"Standard_Parallel_1\",$firstParallel],PARAMETER[\"Standard_Parallel_2\",$secondParallel],PARAMETER[\"Latitude_Of_Origin\",$clat],UNIT[\"Meter\",1.0]]";
            }
        }   
        else {
         print " This type of projection not yet supported - $outprj not created - continuing\n";
         unlink ("$outprj");
        }
       close OUTPROJ;
                
        if (-e $outprj) {
              print " Projection file generated: $outprj\n\n";
        }

    }


