function cst =  matRad_priocstPrep(cst,PriorityList)
        
	cst(:,6) = cell(1);
    for i = 1:numel(PriorityList.ConstraintList)
        cst{PriorityList.ConstraintList{i}.cstIdx,6}{end+1} = PriorityList.ConstraintList{i}.constraint;
    end

end