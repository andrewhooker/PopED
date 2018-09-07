## Test environments
* macOS 10.13.6 (local), R release version
* ubuntu 14.04 (on Travis CI), R release and devel version
* windows Server 2012 R2 x64 (on AppVeyor), R release version
* win-builder, R release and devel version

## R CMD check results
For macOS, ubuntu and windows server there were no ERRORs, WARNINGs or NOTEs. 

For win-builder release and devel versions there was one note:
   
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Andrew C. Hooker <andrew.hooker@farmbio.uu.se>'

Possibly mis-spelled words in DESCRIPTION:
  Foracchia (24:50)
  Nyberg (23:62)
  al (23:72, 24:63)
  et (23:69, 24:60)
     
The first portion is just stating that I am the maintainer.  

The second portion does not recognize the last names of the first authors
of the papers that describe the methods we implement in the package as well as "et al.". 
   
## Downstream dependencies
I have also run R CMD check on downstream dependencies of PopED 
(https://github.com/andrewhooker/PopED/tree/master/revdep). 
All packages passed.