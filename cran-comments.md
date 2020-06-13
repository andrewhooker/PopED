## Test environments
* macOS 10.15.5 (local), R release version
* Ubuntu 16.04.6 LTS (on Travis CI), R release and devel version
* windows Server 2012 R2 x64 (on AppVeyor), R release version
* win-builder, R release and devel version

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

For win-builder release and devel versions there was one
Note_to_CRAN_maintainers on the win-builder versions:
   
* checking CRAN incoming feasibility ... Note_to_CRAN_maintainers
Maintainer: 'Andrew C. Hooker <andrew.hooker@farmbio.uu.se>'

This is just stating that I am the maintainer.  

## Reason for submission
A number of platforms on CRAN check were giving warnings about

* Missing or unexported objects:
     'dplyr::rbind_all' 'dplyr::rbind_list' 

* Found the following files/directories:
     'PopED_output_summary_D_cont_opt_1.txt' 'PopED_output_summary_RS.txt'
     'PopED_output_summary_mfea_opt_1.txt' 'test.csv' 

* Namespace in Imports field not imported from: ‘tidyr’
     All declared Imports should be used.
     

All of these problems have been fixed in the current version. 
Plus many new features have been added.
   
## Downstream dependencies
I have also run R CMD check on downstream dependencies of PopED (currently: ncappc). The package passed.
