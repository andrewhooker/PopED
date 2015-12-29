## Test environments
* local OS X (10.11.2) install, R 3.2.3
* ubuntu 12.04 (on travis-ci), R 3.2.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

## Downstream dependencies
I have also run R CMD check on downstream dependencies of PopED 
(https://github.com/wch/checkresults/blob/master/httr/r-release). 
All packages that I could install passed except:
  
  * Ecoengine: this appears to be a failure related to config on 
that machine. I couldn't reproduce it locally, and it doesn't 
seem to be related to changes in httr (the same problem exists 
                                       with httr 0.4).