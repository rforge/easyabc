
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<h1>EasyABC: a R package to perform efficient approximate Bayesian computation sampling schemes</h1>

<h2>About EasyABC</h2>

<p>The package EasyABC enables to launch a series of simulations of a computer code from the R platform, and to retrieve the simulation outputs in an appropriate format for post-processing treatments.</p>

<p>Four sequential sampling schemes and three coupled-to-MCMC schemes are implemented. EasyABC further enables to launch the simulations in parallel on multiple cores of a computer.</p>

<p><b>Package maintainer:</b> <a href="mailto:nicolas.dumoulin@irstea.fr">Nicolas Dumoulin</a></p>
<p><b>Developers:</b> <a href="mailto:franck.jabot@irstea.fr">Franck Jabot, <a href="mailto:thierry.faure@irstea.fr">Thierry Faure</a>, <a href="mailto:nicolas.dumoulin@irstea.fr">Nicolas Dumoulin</a></p>

<p>Development of EasyABC has been supported by the Irstea project DynIndic and by the French National Research Agency (ANR) within the SYSCOMM project DISCO (ANR-09-SYSC-003).</p>

<p>Thanks to R-Forge for hosting the project.</p>

<h2>News</h2>
<ul>
  <!--li>11/13/2012 EasyABC 1.0 is now available on CRAN.</li-->
  <li>11/13/2012 EasyABC is now hosted on <a href="http://r-forge.r-project.org/projects/easyabc/">R-Forge</a>.</li>
</ul>

<h2>Current features</h2>
<ul>
  <li>Launching of computer simulations to perform sequential or coupled-to-MCMC sampling schemes for Approximate Bayesian Computation.</li>
  <li>EasyABC currently implements four sequential algorithms:
  <ul>
    <li>Beaumont et al. Biometrika (2009)</li>
    <li>Drovandi &amp; Pettit Biometrics (2011)</li>
    <li>Del Moral et al. Statistics and Computing (2012)</li>
    <li>Lenormand et al. ArXiv (2012)</li>
  </ul></li>
  <ul>
    <li>EasyABC also implements three coupled-to-MCMC schemes:
    <li>Marjoram et al. Pnas (2003)</li>
    <li>Wegmann et al. Genetics (2009)</li>
    <li>A mix of Marjoram and Wegmann's algorithms (Marjoram with calibration step – Wegmann without PLS)</li>
  </ul></li>
  <li>EasyABC additionally enables to launch computer simulations in parallel on multiple cores of a multi-core computer.</li>
  <li>EasyABC also contains two R wrappers for binary codes to ease their launching from the R platform.
</ul></li>

<h2>Obtaining EasyABC</h2>

<p>Version 1.0 is available on CRAN. Simply type install.packages("EasyABC") from within R.</p>

<h2>Obtaining help and more information about EasyABC</h2>

<p>There is a package vignette with more information about EasyABC. Simply type vignette("EasyABC") to view within R.</p>

<p>If you have problems or questions about the code, please read the function documentation ( help(EasyABC) ).</p>

<h2>Release history</h2>

<!--p>Version 1.0 released on 13-12-2012.</p-->

</body>
</html>
