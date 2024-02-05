classdef CPPS_Element_2d2el < RC_Element_2d1el

% This is the child class of RC_Element_2d1el.m.
% All the protected & public data properties of parent class will be ...
% inherited to this class 


% Element class for 2nd order analysis of a 2-dimensional structure
    properties (Access = protected)
        kg_local
        k_local
        f_global

        %% Property Definitions 
        %  Note: Variables defined in ud_2d2el.m are not repeated here
        %
        %  kg_local  == Element's local geometric stiffness matrix            
        %
        %  k_local == Element's local stiffness matrix which is the
        %             addition of the elastic and geometric stiffness
        %             matrix
        %  f_global == The global forces for each member which were
        %              recovered locally using the natural deformation
        %              approach and subsequently trasnformed to global
        %              coordinates 

        
    end

  methods (Access = public)
      % Constructor
        function self = CPPS_Element_2d2el(element_nodes, A, Ayy, Izz, E,...
                v, truss)

            self = self@RC_Element_2d1el(element_nodes, A, Ayy, Izz, E, v,...
                truss);
            self.f_local = zeros(6,1); % Initialize the force vector
            self.ComputeElementGeometricStiffMatrix();
        end

        %% GetFGlobal
        % Retrieves the elements global internal forces for access in the
        % analysis class error calculation
        function f_global = GetFGlobal(self)
            f_global = self.f_global;
        end  
        
       %% Compute Forces 
       % This function compuetes the forces for each element in global
       % coordiantes. The forces are recovered using the natural
       % deformation approach. 
        function ComputeForces(self, del_global)

            % Call the del_global and del_local object created in the 
            % RC_Element_2d1el class
            self.del_global = del_global;
            self.del_local = self.gamma * self.del_global; % Multiply by...
            % gamma since del_local is being updated at every step
            
            % Apply Natrual Transformation Approach 
            un = (self.del_local(4)-self.del_local(1)) +...
                ((self.del_local(5)-self.del_local(2))^2 +...
                (self.del_local(4)-self.del_local(1))^2)/(2*self.L);

            thetaR = atan((self.del_local(5)-self.del_local(2))/...
                (self.L + self.del_local(4)-self.del_local(1)));

            thetaAN = self.del_local(3)-thetaR;
            thetaBN = self.del_local(6)-thetaR;
            deltaN = [0;0;thetaAN;un;0;thetaBN];
            dF = (self.ke_local + self.kg_local)*deltaN;
            self.f_local= self.f_local + dF;

        end
        
        %% Update Transformation Matrix
        % This function updates the transformation matrix of the element
        % for each step.
        function UpdateTransformationMatrix(self)
            self.ComputeTransformationMatrix()
        end

        %% Compute Element Geometric Stiffness Matrix
        % Computes the geoemetric stiffness matrix for the element
        function ComputeElementGeometricStiffMatrix(self)
            
            L = self.L;
            P = self.f_local(4);
  
            kg = (P/L) * [1    0      0    -1    0       0;...
                          0   6/5   L/10    0  -6/5     L/10;...
                          0  L/10  2*L^2/15 0  -L/10  -L^2/30;...
                         -1    0      0     1    0       0;...
                          0  -6/5  -L/10    0   6/5    -L/10;...
                          0  L/10  -L^2/30  0  -L/10  2*L^2/15];
            
            self.kg_local = sparse(kg);      
        end

        % Compute Global Stiffness        
        function ComputeGlobalStiffnessMatrix(self)
            self.k_local = self.ke_local + self.kg_local;
            self.k_global = self.gamma'*self.k_local*self.gamma;  
        end
        
        %% Compute Global Element Forces
        % Computes the elements global forces based on updated geometry
        function ComputeGlobalElementForces(self)
            self.f_global = self.gamma'*self.f_local;
        end   
    end    
end
