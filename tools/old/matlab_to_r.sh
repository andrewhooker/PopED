#!/bin/csh
cp $1 $2
ex -s $2 <<eof
   g/%/s//#/g
   g/function\(..*\)=\(..*\)(\(..*\)/s//\2 <-function( \3 { \1/
   g/end/s//   } #/
   g/for\(..*\)=\(..*\):\(..*\)/s//for ( \1 in \2 : \3 ) {/
   g/_/s//./g
   g/;/s///g
   g/==/s//@@/g
   g/=/s//<-/g
   g/@@/s//==/g
   g/zeros(/s//matrix(0,/g
   g/ones(/s//matrix(1,/g
   g/eye(/s//diag(1,/g
   g/\/s//solve(,)/g
   g/fsolve('\(..*\)'/s//ms(~\1 /g
   g/param(\(..*\))/s//param[ \1 ] /g
   g/var(\(..*\))/s//var[ \1 ] /g
   g/mod1(\(..*\)/s//mod1[ \1 /g
   wq
eof