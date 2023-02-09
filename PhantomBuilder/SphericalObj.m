classdef SphericalObj < VolumeObj
    properties
        radius;
    end
    methods
        function obj = SphericalObj(name,type,DoseObjectiveFunction,offset,radius)
            obj@VolumeObj(name,type,DoseObjectiveFunction,offet);
            obj.radius = radius;
        end
        
        function [cst] = initializeParameters(obj,ct,cst)
            cst = initializeParameters@VolumeObj(obj,cst);
            center = round(ct.cubeDim/2);
            VOIHelper = zeros(ct.cubeDim);
            offsets = obj.offset;

            for x = 1:ct.cubeDim(1)
                for y = 1:ct.cubeDim(2)
                   for z = 1:ct.cubeDim(3)
                      currPost = [x y z]  + offsets - center;
                      if  (sqrt(sum(currPost.^2)) < obj.radius)
                            VOIHelper(y,x,z) = 1;
                      end
                   end
                end
                x
            end
            
            cst{obj.idx,4}{1} = find(VOIHelper);
            
        end
    end
end