classdef matRad_OptimizerGamultiobj < matRad_Optimizer
% matRad_OptimizerFmincon implements the interface for the fmincon optimizer 
% of the MATLAB Optiization toolbox
%    
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties
        options     %the optimoptions for gamulitobj
        xs          %last optimization result
        resultInfo  %info struct about last results
        wResult
    end
    
    methods
        function obj = matRad_OptimizerGamultiobj
            %matRad_OptimizerGamultiobj
            %   Construct an instance of the gamultiobj optimizer from the Optimization Toolbox           
            matRad_cfg = MatRad_Config.instance();
            obj.wResult = [];
            obj.xs = [];
            obj.resultInfo = [];
            
            %createDefaultOptimizerOptions Constructs a set of default
            %options for the optimizer to use
            obj.options = optimoptions('gamultiobj',...
                'FunctionTolerance',1e-4,...
                'Display', 'diagnose',...
                'PlotFcn','gaplotpareto',...
                'UseParallel',true);
            end
                
        function obj = optimize(obj,w0,optiProb,dij,cst)
            %optimize carries Out the optimization
            % obtain lower and upper variable bounds
            lb = optiProb.lowerBounds(w0(1,:));
            ub = optiProb.upperBounds(w0(1,:));
            obj.options.InitialPopulationMatrix = w0';            
            % Informing user to press q to terminate optimization
            %fprintf('\nOptimzation initiating...\n');
            %fprintf('Press q to terminate the optimization...\n');
            
            % Run gamultiobj.
            numel(w0)
            [obj.xs,fVal,exitflag,info] = gamultiobj(@(x) obj.gamultiobj_objAndGradWrapper(x,optiProb,dij,cst),...
                size(w0,1),... % Number of variables
                [],[],... % Linear Constraints we do not explicitly use
                [],[],... % Also no linear inequality constraints
                lb,ub,... % Lower and upper bounds for optimization variable
                @(x) obj.gamultiobj_nonlconWrapper(x,optiProb,dij,cst),...
                obj.options); % Non linear constraint structure);
            
            obj.resultInfo = info;
            obj.resultInfo.fVal = fVal;
            obj.resultInfo.exitflag = exitflag;
        end
        
        function fs = gamultiobj_objAndGradWrapper(obj,x,optiProb,dij,cst)
            x = transpose(x);
            fs = optiProb.matRad_objectiveFunctions(x,dij,cst);
        end
        
        function [c,cEq] = gamultiobj_nonlconWrapper(obj,x,optiProb,dij,cst)
            x = transpose(x);
            %Get the bounds of the constraint
            [cl,cu] = optiProb.matRad_getConstraintBounds(cst);
                    
            %Get finite bounds
            clFinIx = isfinite(cl);
            cuFinIx = isfinite(cu);
            
            % Some checks
            assert(isequal(size(cl),size(cu)));
            assert(all(cl <= cu));
            
            %For fmincon we need to separate into equalty and inequality
            %constraints
            isEqConstr = (cl == cu);
            eqIx = isEqConstr;
            ineqIx = ~isEqConstr;
            
            %Obtain all constraint functions and derivatives
            cVals = optiProb.matRad_constraintFunctions(x,dij,cst);   
            
            %Subselection of equality constraints
            cEq = cVals(eqIx & clFinIx); %We can only rely on cl indices here due to the equality index
                        
            %Prepare inequality constraints:
            %We need to separate upper and lower bound constraints for
            %gamultiobj
            cL = cl(ineqIx & clFinIx) - cVals(ineqIx & clFinIx);
            cU = cVals(ineqIx & cuFinIx) - cu(ineqIx & cuFinIx);
            
            %build the inequality jacobian
            c = [cL; cU];
        end
        
        function [statusmsg,statusflag] = GetStatus(obj)
            try 
                statusmsg = obj.resultInfo.message;
                if obj.resultInfo.exitflag == 0
                    statusflag = 0;
                elseif obj.resultInfo.exitflag > 0
                    statusflag = 1;
                else 
                    statusflag = -1;
                end
            catch
                statusmsg = 'No Last Optimizer Status Available!';
                statusflag = -1;
            end
        end
    end
    
    methods (Static)    
        function available = IsAvailable()
            %'gamultiobj' is a p-code file in the optimization toolbox
            available = exist('gamultiobj') == 6;
        end
    end
end