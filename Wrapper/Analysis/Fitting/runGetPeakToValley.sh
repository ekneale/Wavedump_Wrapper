( sleep 10 ; echo .q  ) | root -l 'GetPeakToValley.C++("./outputFile.root",1)';
open ./Plots/*.pdf;