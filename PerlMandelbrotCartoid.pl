#! /usr/bin/perl
#
# This little perl script computes the locations of the 
# buds along the perimiter of the mandelbrot cardiod.
# 

for ($i=2; $i<=12; $i++) 
{

   	$t = 1/$i;
	$p = 2*3.14159265358979*$t;

	$x = 0.5 * cos ($p) - 0.25 * cos (2*$p);
	$y = 0.5 * sin ($p) - 0.25 * sin (2*$p);

	printf  "bud %3d is at theta=%8.6f   re(z)=%8.6f   im(z)=%8.6f\n", $i, $p, $x, $y;
}
