function [cube, metadata] = matRad_readCube(filename)
% matRad Cube read wrapper
% determines the extension and assigns the appropriate reader to it
% 
% call
%   matRad_readCube(filename)
%
% input
%   filename:   full path of the file
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Instantiate matrad config
matRad_cfg = MatRad_Config.instance();

if ~exist(filename,'file')
    error(['File ' filename ' does not exist!']);
end

[pathstr,name,ext] = fileparts(filename);
if strcmp(ext,'.gz')
    [~,name,subext] = fileparts(name);
    ext = strcat(subext,ext);
end

switch ext
    case {'.nrrd','.NRRD'}
        matRad_cfg.dispInfo('Reading NRRD: "%s" ...',filename);
        [cube, metadata] = matRad_readNRRD(filename);
        matRad_cfg.dispInfo('Done!\n');
    case {'.nii','.nii.gz'}
        matRad_cfg.dispInfo('Reading NifTI: "%s" ...',filename);
        [cube, metadata] = matRad_readNifTI(filename);
        matRad_cfg.dispInfo('Done!\n');
    otherwise
        matRad_cfg.dispError('Extension %s not (yet) supported!',ext);
end
metadata.name = name;
metadata.path = pathstr;m

end

