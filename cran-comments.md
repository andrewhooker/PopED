## Changes

Changes in this version of PopED are:

* A fix to a bug pointed out by Duncan Murdoch about the use of 
  C-like `return` code.
  
* A change to the PopED NAMESPACE so that an S3 method can be seen by users.

## Test environments
* local macOS (10.12.1) install, R 3.3.2
* ubuntu 12.04 (on travis-ci), R 3.3.1
* win-builder (release and devel)

## R CMD check results
For macOS and ubuntu there were no ERRORs, WARNINGs or NOTEs. 

For win-builder release and devel versions there was one note:
   
* checking CRAN incoming feasibility ... NOTE
    + Maintainer: 'Andrew C. Hooker <andrew.hooker@farmbio.uu.se>'
    + Possibly mis-spelled words in DESCRIPTION:
      pharmacometric (19:9)
     
The first portion is just stating that I am the maintainer.  

The second portion does not recognize the well known word "pharmacometric",
see for example: https://en.wikipedia.org/wiki/Pharmacometrics
   
## Downstream dependencies
There are currently no downstream dependencies for this package.
