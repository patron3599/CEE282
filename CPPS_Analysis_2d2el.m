classdef CPPS_Analysis_2d2el < RC_Analysis_2d1el

% This is the child class of RC_Analysis_2d1el.m.
% All the protected & public data properties of parent class will be ...
% inherited to this class 


% Analysis class for 2nd order analysis of a 2-dimensional structure
    

    properties (Access = protected)
        maxsteps      
        numsteps
        restart
        apratios
        limit_state
        h_stat_mes
        dconcen
        ratio_req
        stop_ratio
        dconcen_t
        step
        DEFL_STEP
        REACT_STEP
        APRATIOS
        LIMIT_STATE
        R
        E
        load_norm_E
        energy_norm_E
        
        %% Property Definitions 
        %  Note: Variables defined in ud_2d2el.m are not repeated here
        %
        %  maxsteps  == Is the maximum number of steps that the code will ...
        %               run. This is determined by selecting the minumum ...
        %               of numsteps or the ratio of stop_ratio/ratio_req        
        %
        %  dconcen == The incremental concentraed loads which are stored
        %              along with the full concentrated loads
        %
        %  step == An index to indicate the current step in which the
        %          analysis is in 
        %
        %  R  == A vector with dimensions (ndof,1) which is used to store
        %        the calculated reactions after each step 
        %
        %  E == A vector with dimensions (freedof,1) which is used to store
        %       the calcualted error for each step 
        %
        %  load_norm_E, energy_norm_E == Vectors with dimensions (step,1)
        %                                which are used to the energy ...
        %                                indices for the current step
        %  DEFL_STEP  == Similar to DEFL, however this property only stores
        %                values for the current step and is updated with...
        %                each subsequent one 
        %
        %  REACT_STEP  == Similar to REACT, however this property only...
        %                stores values for the current step and is updated...
        %                with each subsequent one 
        %

    end
     methods (Access = public)
        function self = CPPS_Analysis_2d2el(nnodes, coord, fixity, concen,...
                nele, ends, A, Ayy, Izz, E, v, truss,numsteps,ratio_req,...
                stop_ratio,restart,apratios,limit_state,h_stat_mes)
            
            self = self@RC_Analysis_2d1el(nnodes, coord, fixity, concen,...
                nele, ends, A, Ayy, Izz, E, v, truss);

            self.ratio_req=ratio_req;
            self.stop_ratio=stop_ratio;
            self.restart=restart;
            self.dconcen = concen.*ratio_req;
            self.dconcen_t = self.dconcen';
            self.apratios=apratios;
            self.limit_state=limit_state;
            self.h_stat_mes=h_stat_mes;
            self.numsteps=numsteps;
            self.maxsteps=min(numsteps,floor(stop_ratio/ratio_req));
            
        end

        %% Run Analysis
        % Runs the 2nd order analysis
        function RunAnalysis(self)
            self.InitializeOutputVariables()
            self.CreateStiffnessMatrix();
            self.CreateLoadVectors();
            self.step = 1;
            
            % Runs the analysis as long as AFLAG is not activated 
            if self.AFLAG
                while self.step <= self.maxsteps && self.LIMIT_STATE == 0

                    self.APRATIOS = [self.APRATIOS;self.step*...
                        self.ratio_req];

                    self.runIteration()
                    self.step = self.step + 1;
                end
            end
            
            RC_Plot_Errors(self.load_norm_E,self.energy_norm_E,...
                self.APRATIOS);
            
        end
        
        %% Run Iteration
        % Runs each iteration of the 2nd order analysis
        function runIteration(self)

            self.ComputeDisplacementsReactions(self.step);
            self.RecoverElementForces(self.step) 
            self.ReviseNodes();
            self.UpdateStiffnessMatrices();
            self.CreateStiffnessMatrix();
            self.CheckLimitState();
            self.ComputeError();
        end
        
        %% Get Mastan2 Returns
        %  Returns the matrices that need to be returned to Mastan2
        function [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE] =...
                GetMastan2Returns(self)

            DEFL = self.DEFL;
            REACT = self.REACT;
            ELE_FOR = self.ELE_FOR;
            AFLAG = self.AFLAG;
            APRATIOS = self.APRATIOS;
            LIMIT_STATE = self.LIMIT_STATE; 
        end
       
    end
    
    methods (Access = protected)

        %% Initialize Output Variables
        function InitializeOutputVariables(self)
        % Generates the required output variables and also adjusts property
        % parameters such that they can be used for linear indexing.
            % Re-arrange the dimensions of DEFL_STEP and REACT_STEP such
            % that they can be concatenated to DEFL and REACT, respectively
            self.DEFL_STEP = zeros(self.num_dof_node, self.nnodes);
            self.REACT_STEP = zeros(self.num_dof_node, self.nnodes);
            
            % Initialize for first step of analysis 
            self.DEFL = zeros(self.nnodes,self.num_dof_node);
            self.REACT = zeros(self.nnodes,self.num_dof_node);
            self.ELE_FOR = zeros(self.nele,self.num_dof_node*2);
            self.APRATIOS = self.apratios; 
            self.E = []; 
            self.load_norm_E = [];
            self.energy_norm_E = [];
            self.LIMIT_STATE = 0; 
           
        end

        %% Create Nodes
        % Create the (nnodes, 1) vector of node objects representing all .. 
        % the nodes in the structure as a second order object. This method...
        % was edited from the RC_Analysis_2D1el class.
        function CreateNodes(self)
            for i = 1:self.nnodes             
                self.nodes = [self.nodes;...
                    CPPS_Node_2d2el(i, self.coord_t(:,i))];
            end
        end

        %% Revise Nodes
        % Changes the node's coordinate based off the deflections computed
        % in the previous step
        function ReviseNodes(self)
            for i = 1:self.nnodes
                DEFL_STEP_t = self.DEFL_STEP';
                self.nodes(i).UpdateNodeCoord(DEFL_STEP_t(i,[1,2]))
            end
        end
        
        %% Create Elements
        %  Creates the (nele, 1) vector of element objects representing all...
        % the elements in the structure a 2nd order object. This method was
        % edited from the RC_Analysis_2D1el class
        function CreateElements(self, A, Ayy, Izz, E, v)

            for i = 1:self.nele
                self.elements = [self.elements; ...
                    CPPS_Element_2d2el(self.nodes(self.ends(i, 1:2)),...
                    A(i), Ayy(i), Izz(i), E(i), v(i), self.truss)];
            end
        end
        
        %% Create Load Vectors
        %  Updates the applied load vectors for each iteration. This method
        %  was edited from the RC_Analysis_2D1el class
        function CreateLoadVectors(self)
            
            % Compute vector of concentrated loads applied at the free and...
            % support degrees of freedom using ...
            % linear indexing of the "concen_t" matrix
            self.Pf = self.dconcen_t(self.dof_free);
            self.Psupp = self.dconcen_t(self.dof_supp);
            
            % Compute vector of specified displacements using...
            % linear indexing of the "fixity_t" matrix
            self.deln = self.fixity_t(self.dof_disp);
            
        end
        
        %% Compute Displacements Reactions
        %  Compute the displacements and reactions and format them to
        %  return to Mastan2. This method was edited from the...
        % RC_Analysis_2D1el class.
        function ComputeDisplacementsReactions(self,i)
            
            % Compute the displacements
            self.delf = self.Kff \ (self.Pf - self.Kfn*self.deln);
            
            % Compute the reactions, accounting for loads applied ...
            % directly on the supports
            self.Ps = self.Ksf*self.delf + self.Ksn*self.deln - self.Psupp;
            self.Pn = self.Knf*self.delf + self.Knn*self.deln;
            
            % Format the computed displacements using linear...
            % indexing of the "DEFL" matrix
            self.DEFL_STEP(self.dof_free) = self.delf;
            self.DEFL_STEP(self.dof_disp) = self.deln;
            if i == 1
                self.DEFL(:,:,i) = self.DEFL_STEP';
            else
                self.DEFL = cat(3,self.DEFL,self.DEFL(:,:,i-1)+...
                    self.DEFL_STEP');
            end
            
            
            % Format the computed reactions using linear indexing of ...
            % the "REACT" matrix
            self.REACT_STEP(self.dof_supp) = self.Ps;
            self.REACT_STEP(self.dof_disp) = self.Pn;
            if i == 1
                self.REACT(:,:,i) = self.REACT_STEP';
            else
                self.REACT = cat(3,self.REACT,self.REACT(:,:,i-1)+...
                    self.REACT_STEP');
            end    
        end
        
        %% Recover Element Forces
        %  Recovers the local element forces of each iteration and formats...
        % them to return to Mastan2. This method was edited from the...
        % RC_Analysis_2D1el class.
        function RecoverElementForces(self, j)
            % Initialize the ELE_FOR matrix with zeros for the first step
            if j ~= 1
            self.ELE_FOR = cat(3,self.ELE_FOR,zeros(self.nele,...
                self.num_dof_node*2));
            end

            for i = 1:self.nele
                
                % Obtain the displacements at the degrees of freedom ...
                % corresponding to element i using linear indexing of the...
                % "DEFL_STEP" matrix
                if j==1
                self.elements(i).ComputeForces(...
                    self.DEFL_STEP(self.elements(i).GetElementDOF()));
                else
                self.elements(i).ComputeForces(...
                    self.DEFL_STEP(self.elements(i).GetElementDOF()));    
                end
                self.ELE_FOR(i,:,j) = self.elements(i).GetFLocal();
            end
        end
        
        
        %% Update Stiffness Matrices
        % Updates the transformation and geometric stiffness matricies to
        % be used in the next iteraton.
        function UpdateStiffnessMatrices(self)
            for i = 1:self.nele
                self.elements(i).UpdateTransformationMatrix();
                self.elements(i).ComputeElementGeometricStiffMatrix();
            end
        end
        
        %% Compute Error
        % Computes the error for all free degrees of freedom and stores the
        % value for each in the created R and E vectors
        function ComputeError(self)
            self.R = zeros(self.num_dof_total,1); % Initialize R 

            for i = 1:self.nele
                self.elements(i).ComputeGlobalElementForces(); 
                f_global = self.elements(i).GetFGlobal;
                element_dof = self.elements(i).GetElementDOF();
               
               % Assemble the R vector for every degree of freedom
               self.R(element_dof) = self.R(element_dof) + f_global;
            end

            applied_load = self.APRATIOS(self.step)*self.concen_t(:);
            E_total = applied_load-self.R;
            
            %The error vector for the current step at free dof
            E_step = E_total(self.dof_free); 
            
            % Concatenate to the E vector from the previous step. 
            self.E = horzcat(self.E,E_step);  
                    
            % Load norm error index for current step 
            applied_load_free = applied_load(self.dof_free);
            load_norm_E = norm(E_step)/norm(applied_load_free);
            self.load_norm_E = [self.load_norm_E;load_norm_E];
            
            % Energy norm error index for current step
            energy_norm_E = (abs(E_step)'*abs(self.delf))/...
                (abs(applied_load_free)'*abs(self.delf));

            self.energy_norm_E = [self.energy_norm_E;energy_norm_E];
        end
            
        %% Check Limit State
        % Checks the limit state of the current Kff matrix 
        function CheckLimitState(self)
            [~, p] = chol(self.Kff);
            if p ==0
                self.LIMIT_STATE = 0; % Structure is loading 
            else
                self.LIMIT_STATE = 1; % Structure is unloading 
            end
        end
    end
end
