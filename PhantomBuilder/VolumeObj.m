classdef (Abstract) VolumeObj
    properties
        idx
        name;
        type;
        TissueClass = 1;
        alphaX = 0.1000;
        betaX = 0.0500;
        Priority = 1;
        Visible = 1;
        visibleColor = [0 0 0];
        center = [0,0,0]; %center of objective
        DoseObjectiveFunction;
        
    end
    methods (Static, Access = private)
        function oldValue = getOrIncrementCount(increment)
        % Private function to manage the counter
            persistent VALUE
            if isempty(VALUE)
                VALUE = 0;
            end
            oldValue = VALUE;
            if nargin > 0
                VALUE = VALUE + increment;
            end
        end 
    end
     methods (Static)
        function value = getInstanceCount()
        % Public access to the counter cannot increment it
            value = cldef.getOrIncrementCount();
        end
    end
    methods
        function obj = VolumeObj(name,type)
        % Increment the counter in the constructor
            VolumeObj.getOrIncrementCount(1);
            obj.idx = VolumeObj.getOrIncrementCount();
            obj.name = name;
            obj.type = type;
        end
        
        function cst = initializeParameters(obj,cst)
            %initialize cst file
            cst{obj.idx,1}                = obj.idx-1;
            cst{obj.idx,2}                = obj.name;
            cst{obj.idx,3}                = obj.type;
            cst{obj.idx,5}.TissueClass    = obj.TissueClass;
            cst{obj.idx,5}.alphaX         = obj.alphaX;
            cst{obj.idx,5}.betaX          = obj.betaX;
            cst{obj.idx,5}.Priority       = obj.Priority;
            cst{obj.idx,5}.Visible        = obj.Visible;
            cst{obj.idx,5}.visibleColor   = obj.visibleColor;
                            
        end
    end
end