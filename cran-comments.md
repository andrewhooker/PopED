## Test environments
* macOS 11.3.1 (local), R-release
* Github actions, windows-latest,  R-release
* Github actions, macOS-latest,  R-release 
* Github actions, ubuntu-20.04,  R-release and R-devel 
* rhub, Debian Linux,  R-devel
* rhub, Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results
There were no ERRORs, WARNINGs.

On some flavors I get the NOTE:

  * checking CRAN incoming feasibility ... 
    Maintainer: 'Andrew C. Hooker <andrew.hooker@farmaci.uu.se>'

Which, of course, is just my e-mail address. Note that this has changed
after a move to a different department at my university.

## Reason for submission
I got a message from Brian Ripley that the PopED-Ex.Rout document in the 
latest development candidate of R
has the following warning:

  Warning in matrix( ... ... ... )
   data length differs from size of matrix
 
There was a bug in a function and has been fixed in the current version. 

A number of new features have also been added.

NOTE: I have changed my e-mail address after a move to a different department 
at my university.
   
## Downstream dependencies
I have also run R CMD check on downstream dependencies of 
PopED (currently: ncappc). The package was not affected by the changes.
