## Test environments
* macOS 11.3.1 (local), R-release (4.1.0)
* Github actions, windows-latest, R-release (4.1.0)
* Github actions, macOS-latest, R-release (4.1.0)
* Github actions, ubuntu-20.04, R-release (4.1.0)
* Github actions, ubuntu-20.04, R-devel (2021-05-19 r80339) 
* rhub, Ubuntu Linux 20.04.1 LTS, R-release (4.1.0), GCC
* rhub, Fedora Linux, R-devel, clang, gfortran
* winbuilder, R-devel (2021-05-18 r80323)

## R CMD check results
There were no ERRORs, WARNINGs.

* From rhub and winbuilder:

checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Andrew C. Hooker <andrew.hooker@farmaci.uu.se>’

New maintainer:
  Andrew C. Hooker <andrew.hooker@farmaci.uu.se>
Old maintainer(s):
  Andrew C. Hooker <andrew.hooker@farmbio.uu.se>

I have changed my e-mail address after a move to a different department 
at my university.

## Reason for submission

* This is a re-submission. 
On the original submission I had:

1 NOTE in the Windows and Debian check of the development version relating to 
DOIs and URLs that have been fixed.

1 WARNING in the Debian devel version about vignette warnings related to 
right_join without the by argument and a lack of library(ggplot2) where needed.
These have been fixed.

* Reason for the original submission:

I got a message from Brian Ripley that the PopED-Ex.Rout document in the 
latest development candidate of R
has the following warning:

  Warning in matrix( ... ... ... )
   data length differs from size of matrix
 
There was a bug in a function and has been fixed in the current version. 

A number of new features have also been added.

* NOTE: I have changed my e-mail address after a move to a different department 
at my university.
   
## Downstream dependencies
I have also run R CMD check on downstream dependencies of 
PopED (currently: ncappc). The package was not affected by the changes.
