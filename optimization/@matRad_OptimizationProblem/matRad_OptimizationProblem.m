classdef matRad_OptimizationProblem < handle
    %matRad_OptimizationProblem Main Class for fluence optimization problems
    % Describes a standard fluence optimization problem by providing the 
    % implementation of the objective & constraint function/gradient wrappers
    % and managing the mapping and backprojection of the respective dose-
    % related quantity
    %
    % References
    %   [1] https://doi.org/10.1093/imanum/draa038
    %   [2] https://doi.org/10.1002/mp.14148
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
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
        BP
        normalizationScheme = struct('type','none');
        objectives = {}; %cell array storing all objectives, has to be initialized at the start
        constraints = {}; %
        objidx;
        constridx;
        quantityOpt = '';
        useMaxApprox = 'logsumexp'; %'pnorm'; %'logsumexp'; %'none';
        p = 30; %Can be chosen larger (closer to maximum) or smaller (closer to mean). Only tested 20 >= p >= 1

        minimumW = NaN;
        maximumW = NaN;
    end

    properties ()

    end
    
    methods
        function obj = matRad_OptimizationProblem(backProjection,cst)
            
            obj.BP = backProjection; %needs to be initalized to have access to setBiologicalDosePrescriptions
            if nargin == 2
                obj.matRad_extractObjectivesFromcst(cst);
                obj.matRad_extractConstraintsFromcst(cst);
            end
        end       
        
        %Objective function declaration
        fVal = matRad_objectiveFunction2(optiProb,w,dij,cst)   
        
        %Objective gradient declaration
        fGrad = matRad_objectiveGradient2(optiProb,w,dij,cst)
        
        %Constraint function declaration
        cVal = matRad_constraintFunctions(optiProb,w,dij,cst)
        
        %Constraint Jacobian declaration
        cJacob = matRad_constraintJacobian(optiProb,w,dij,cst)
        
        %Jacobian Structure
        jacobStruct = matRad_getJacobianStructure(optiProb,w,dij,cst)
        
        [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
        
        function lb = lowerBounds(optiProb,w)
            minW = optiProb.minimumW;
            if isnan(minW)
                lb = zeros(size(w));
            elseif isscalar(minW)
                lb = minW*ones(size(w));
            elseif isvector(minW) && all(size(minW) == size(w))
                lb = minW;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Minimum Bounds for Optimization Problem could not be set!');
            end
        end
        
        function ub = upperBounds(optiProb,w)
            maxW = optiProb.maximumW;
            if isnan(maxW)
                ub = Inf(size(w));
            elseif isscalar(maxW)
                ub = maxW*ones(size(w));
            elseif isvector(maxW) && all(size(maxW) == size(w))
                ub = maxW;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Maximum Bounds for Optimization Problem could not be set!');
            end
        end

        function normalizedfVals = normalizeObjectives(optiProb,fVals)
            switch optiProb.normalizationScheme.type
                case 'none'
                    normalizedfVals = fVals;
                case 'UL'
                    %maybe check that U and L are defined
                    normalizedfVals = (fVals - optiProb.normalizationScheme.L)./(optiProb.normalizationScheme.U-optiProb.normalizationScheme.L); %might have to check that U and L work!
                otherwise
                    matRad_cfg.dispError('Normalization scheme not known!');
            end
        end

        function normalizedGradient = normalizeGradients(optiProb,Gradient)
            switch optiProb.normalizationScheme.type
                case 'none'
                    normalizedGradient = Gradient;
                case 'UL'
                    %maybe check that U and L are defined
                    normalizedfVals = Gradient./(optiProb.normalizationScheme.U-optiProb.normalizationScheme.L); %might have to check that U and L work!
                otherwise
                    matRad_cfg.dispError('Normalization scheme not known!');
            end
        end

        function updatePenalties(obj,newPen) %does it handle grouping?
            if numel(obj.objectives ~= numel(newPen))
                matRad_cfg.dispError('Number of objectives in optimization Problem not equal to number of new penalties to be set!');
            end
            for i=1:numel(newPen)
                obj.objectives{i,1} = newPen(i);
            end
        end

    end
    
    methods (Access = protected)
        function [val,grad] = logSumExp(optiProb,fVals)
            % [1] stable log sum exp trick
            [fMax,ixMax] = max(fVals(:));
            
            ix = true(numel(fVals),1);
            ix(ixMax) = 0;

            tmp = exp(fVals - fMax);
                       
            expSum = sum(tmp(ix));
            val = fMax + log1p(expSum); %log1p(x) = Matlab's numerically accurate log(1+x) 
            
            grad = tmp ./ (1 + expSum);
        end
        
        function [val,grad] = pNorm(optiProb,fVals,n)
            % Implemented as proposed in [2] including a normalization for stability of the exponent.
            if nargout < 3
                n = numel(fVals);
            end
            
            p = optiProb.p;
            
            valMax = max(fVals(:));            
            tmp = fVals./valMax;            
            
            pNormVal = sum(tmp(:).^p)^(1/p);
            
            fac = (1/n)^(1/p);
            val = valMax*fac*pNormVal;

            grad = fac * (tmp ./ pNormVal).^(p-1);
        end
    end     
    
    
    methods (Access = private)
        function matRad_extractObjectivesFromcst(optiProb,cst)
            %used to store objectives in cell array as property of optimization Problem
            optiProb.objidx = [];
            optiProb.objectives = {};
            
            for i = 1:size(cst,1) % loop over cst
                for j = 1:numel(cst{i,6})
                    %check whether dose objective or constraint
                    objective = cst{i,6}{j};
                    if isstruct(cst{i,6}{j})
                        objective =  matRad_DoseOptimizationFunction.createInstanceFromStruct(objective);
                    end
                    if contains(class(objective),'DoseObjectives')
                        optiProb.objidx = [objidx;i,j];
                            objective = optiProb.BP.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);
                        optiProb.objectives(end+1) = {objective};
                    end
                end
            end
            

        end

        
        function matRad_extractConstraintsFromcst(optiProb,cst) %need to check
            %used to store constraints in cell array as property of optimization Problem 
            optiProb.constridx = [];
            optiProb.constraints = {};
            
            for i = 1:size(cst,1) % loop over cst
                for j = 1:numel(cst{i,6})
                    %check whether dose objective or constraint
                    constraint = cst{i,6}{j};
                    if isstruct(cst{i,6}{j})
                        constraint =  matRad_DoseOptimizationFunction.createInstanceFromStruct(constraint);
                    end
                    if contains(class(constraint),'DoseConstraints')
                        optiProb.constridx = [constridx;i,j];
                        constraint = optiProb.BP.setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);
                        optiProb.constraints(end+1) = {constraint};
                    end
                end
            end
        end
end

