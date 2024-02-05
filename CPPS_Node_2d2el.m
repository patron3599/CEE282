classdef CPPS_Node_2d2el < RC_Node_2d1el

% This is the child class of RC_Node_2d1el.m.
% All the protected & public data properties of parent class will be ...
% inherited to this class 

% Node class for 2nd order analysis of a 2-dimensional structure
        
    methods (Access = public)
        %% Constructor
        function self = CPPS_Node_2d2el(node_num,node_coord)
            self = self@RC_Node_2d1el(node_num,node_coord);
        end

        %% UpdateNodeCoord
        % Updates the node object's location after each step of the
        % analysis by an amount d_delta. d_delta is a 1x2 vector which ...
        % contains the change in displacement for that step 

        function UpdateNodeCoord(self,d_delta)
            self.node_coord = self.node_coord + d_delta';    
        end

        
    end   
end
