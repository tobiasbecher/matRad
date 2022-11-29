classdef SphericalObj < VolumeObj
    properties
        radius
    end
    methods
        function obj = SphericalObj()
        obj@VolumeObj();
        end
    end
end