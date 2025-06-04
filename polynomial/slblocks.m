function blkStruct = slblocks
%SLBLOCKS  Define the block library for a specific Toolbox or Blockset

%	Author(s): M. Hromcik, 6-12-2000
%	Copyright (c) 2000 by Polyx, Ltd.
%  Version 3

blkStruct.Name = ['Polynomial' sprintf('\n') 'Toolbox 3.0'];
blkStruct.OpenFcn = 'polblock';
blkStruct.MaskDisplay = 'disp(''FRAC'')';%'fprintf('' Matrix \n Fraction'')';

%end .. slblocks


