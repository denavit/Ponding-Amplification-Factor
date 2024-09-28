import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, ceil, pi
from PyPonding import FE, PondingLoadCell
from PyPonding import opensees as ops
from libdenavit import camber, OpenWebSteelJoist, JoistGirder
from libdenavit.OpenSees import AnalysisResults
from sji_load_tables import get_joist_data

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

class WarehouseBay:
    """
    Class to perform ponding analyses on an idealized rectangular bay. 
    
    
    """
    
    def __init__(self, **attrs):
    
        # Geometry
        self.L12 = attrs['L12']
        self.LAB = attrs['LAB']
        self.LBC = attrs['LBC']
        self.LCD = attrs['LCD']
        self.number_of_joist_spaces = attrs['number_of_joist_spaces']

        # Loads on bay
        self.dead_load_uniform  = attrs['dead_load_uniform'] 
        self.snow_density       = attrs['snow_density'] 
        self.snow_height        = attrs['snow_height'] 
        self.water_density      = attrs['water_density'] 
        
        # Loads on primary members (force per length)
        self.dead_load_on_primary_member_B = attrs.get('dead_load_on_primary_member_B', 0.0)
        self.dead_load_on_primary_member_C = attrs.get('dead_load_on_primary_member_C', 0.0)
        
        # Load factors and options
        self.alpha               = attrs['alpha']
        self.load_factor_dead    = attrs['load_factor_dead']    # Dead
        self.load_factor_ponding = attrs['load_factor_ponding'] # Impounded Water
        self.load_factor_snow    = attrs['load_factor_snow']    # Snow
        self.consider_snow_and_water_overlap = attrs.get('consider_snow_and_water_overlap', True)

        self.run_factored_analysis_after_ponding_analysis = attrs.get('run_factored_analysis_after_ponding_analysis', False)
        self.additional_load_factor_dead    = attrs.get('additional_load_factor_dead', 1.0)
        self.additional_load_factor_ponding = attrs.get('additional_load_factor_ponding', 1.0)
        self.additional_load_factor_snow    = attrs.get('additional_load_factor_snow', 1.0)

        # Top of roof elevation
        self.zA = attrs['zA']
        self.zB = attrs['zB']
        self.zC = attrs['zC']
        self.zD = attrs['zD']

        # Camber
        self.camber_of_primary_member_B = attrs['camber_of_primary_member_B']
        self.camber_of_primary_member_C = attrs['camber_of_primary_member_C']
        self.camber_of_secondary_members_AB = attrs['camber_of_secondary_members_AB']
        self.camber_of_secondary_members_BC = attrs['camber_of_secondary_members_BC']
        self.camber_of_secondary_members_CD = attrs['camber_of_secondary_members_CD']

        # Material and section properties
        self.modulus_of_elasticity = attrs['modulus_of_elasticity']
        self.cross_sectional_area  = attrs['cross_sectional_area']
        self._moment_of_inertia_AB = attrs.get('_moment_of_inertia_AB',None)
        self._moment_of_inertia_BC = attrs.get('_moment_of_inertia_BC',None)
        self._moment_of_inertia_CD = attrs.get('_moment_of_inertia_CD',None)
        self._moment_of_inertia_B  = attrs.get('_moment_of_inertia_B',None)
        self._moment_of_inertia_C  = attrs.get('_moment_of_inertia_C',None)

        # Analysis options
        self.analsis_engine         = attrs.get('analsis_engine', 'FE')
        self.include_ponding_effect = attrs.get('include_ponding_effect', True)
        self.num_ele_secondary      = attrs.get('num_ele_secondary', 20)
        self.num_subcell_X          = attrs.get('num_subcell_X', 4) # Number of ponding sub-cells along joist direction
        self.num_subcell_Y          = attrs.get('num_subcell_Y', 4) # Number of ponding sub-cells along joist girder direction       
        self.MAX_ITER               = attrs.get('MAX_ITER', 50) # Maximum number of ponding analysis iterations
        self.tol                    = attrs.get('tol', 0.00001) # Ponding analysis tolerance
       
    @property
    def moment_of_inertia_AB(self):
        if self._moment_of_inertia_AB is not None:
            return self._moment_of_inertia_AB
        else:
            raise Exception('moment of inertia from joist object not yet implemented')
            
    @property
    def moment_of_inertia_BC(self):
        if self._moment_of_inertia_BC is not None:
            return self._moment_of_inertia_BC
        else:
            raise Exception('moment of inertia from joist object not yet implemented')

    @property
    def moment_of_inertia_CD(self):
        if self._moment_of_inertia_CD is not None:
            return self._moment_of_inertia_CD
        else:
            raise Exception('moment of inertia from joist object not yet implemented')
            
    @property
    def moment_of_inertia_B(self):
        if self._moment_of_inertia_B is not None:
            return self._moment_of_inertia_B
        else:
            raise Exception('moment of inertia from joist object not yet implemented')
            
    @property
    def moment_of_inertia_C(self):
        if self._moment_of_inertia_C is not None:
            return self._moment_of_inertia_C
        else:
            raise Exception('moment of inertia from joist object not yet implemented')
            
    def Run_Analysis(self,water_level):
           
        # Determine profile based on top of roof elevations and camber
        roof_profile = np.zeros((self.number_of_joist_spaces+1,3*self.num_ele_secondary+1))
        for i in range(self.number_of_joist_spaces+1):
            x = i/self.number_of_joist_spaces           
            zA = self.zA
            zB = self.zB + camber(x*self.L12,self.L12,self.camber_of_primary_member_B)
            zC = self.zC + camber(x*self.L12,self.L12,self.camber_of_primary_member_C)
            zD = self.zD

            for j in range(self.num_ele_secondary+1):
                x = j/self.num_ele_secondary
                roof_profile[i,j] = zA + x*(zB-zA) + \
                    camber(x*self.LAB,self.LAB,self.camber_of_secondary_members_AB)
                roof_profile[i,self.num_ele_secondary+j] = zB + x*(zC-zB) + \
                    camber(x*self.LBC,self.LBC,self.camber_of_secondary_members_BC)
                roof_profile[i,2*self.num_ele_secondary+j] = zC + x*(zD-zC) + \
                    camber(x*self.LCD,self.LCD,self.camber_of_secondary_members_CD)
            
        #print(np.transpose(roof_profile))

        # Define ponding load cells
        ponding_load_cells = [[0] * 3*self.num_ele_secondary for i in range(self.number_of_joist_spaces)]
        for i in range(self.number_of_joist_spaces):
            for j in range(3*self.num_ele_secondary):
                
                iCell = PondingLoadCell.PondingLoadCell3d()
                iCell.id      = '%i,%i' % (i,j)
                
                iCell.xI      = i/self.number_of_joist_spaces*self.L12
                iCell.zI      = roof_profile[i,j+1]
                iCell.xJ      = (i+1)/self.number_of_joist_spaces*self.L12
                iCell.zJ      = roof_profile[i+1,j+1]
                iCell.xK      = (i+1)/self.number_of_joist_spaces*self.L12
                iCell.zK      = roof_profile[i+1,j]
                iCell.xL      = i/self.number_of_joist_spaces*self.L12
                iCell.zL      = roof_profile[i,j]

                if j < self.num_ele_secondary:
                    iCell.yI = -(j+1)/self.num_ele_secondary*self.LAB
                    iCell.yJ = -(j+1)/self.num_ele_secondary*self.LAB
                    iCell.yK = -j/self.num_ele_secondary*self.LAB
                    iCell.yL = -j/self.num_ele_secondary*self.LAB
                elif j < 2*self.num_ele_secondary:
                    iCell.yI = -(j+1)/self.num_ele_secondary*self.LBC + self.LAB
                    iCell.yJ = -(j+1)/self.num_ele_secondary*self.LBC + self.LAB
                    iCell.yK = -j/self.num_ele_secondary*self.LBC + self.LAB
                    iCell.yL = -j/self.num_ele_secondary*self.LBC + self.LAB     
                elif j < 3*self.num_ele_secondary:
                    iCell.yI = -(j+1)/self.num_ele_secondary*self.LCD + self.LBC + self.LAB
                    iCell.yJ = -(j+1)/self.num_ele_secondary*self.LCD + self.LBC + self.LAB
                    iCell.yK = -j/self.num_ele_secondary*self.LCD + self.LBC + self.LAB
                    iCell.yL = -j/self.num_ele_secondary*self.LCD + self.LBC + self.LAB
                else:
                    raise Exception('Should not reach here.')

                iCell.gamma   = self.alpha*self.load_factor_ponding*self.water_density
                iCell.na      = self.num_subcell_X
                iCell.nb      = self.num_subcell_Y

                iCell.gammas  = self.snow_density
                iCell.hs      = self.alpha*self.load_factor_snow*self.snow_height
                
                iCell.return_water_load_only = True

                ponding_load_cells[i][j] = iCell
               
        # Run Ponding Analysis
        ponding_load      = np.zeros((self.number_of_joist_spaces+1,3*self.num_ele_secondary+1))
        ponding_load_last = np.zeros((self.number_of_joist_spaces+1,3*self.num_ele_secondary+1))
        for iteration in range(self.MAX_ITER):
            
            # Reset ponding load arrays
            for i in range(self.number_of_joist_spaces+1):
                for j in range(3*self.num_ele_secondary+1):
                    ponding_load_last[i,j] = ponding_load[i,j]
                    ponding_load[i,j] = 0
                    
            # Compute new ponding load
            for i in range(self.number_of_joist_spaces):
                for j in range(3*self.num_ele_secondary):
                    if iteration == 0:
                        ponding_load_cells[i][j].dzI = 0.0
                        ponding_load_cells[i][j].dzJ = 0.0
                        ponding_load_cells[i][j].dzK = 0.0
                        ponding_load_cells[i][j].dzL = 0.0
                    else:
                        ponding_load_cells[i][j].dzI = analysis_results.total_deflection[i,j+1]
                        ponding_load_cells[i][j].dzJ = analysis_results.total_deflection[i+1,j+1]
                        ponding_load_cells[i][j].dzK = analysis_results.total_deflection[i+1,j]
                        ponding_load_cells[i][j].dzL = analysis_results.total_deflection[i,j]
                        
                    f = ponding_load_cells[i][j].get_load_vector(water_level)
                    ponding_load[i,j+1]   += f[0]
                    ponding_load[i+1,j+1] += f[1]
                    ponding_load[i+1,j]   += f[2]
                    ponding_load[i,j]     += f[3]

            # Perform analysis on secondary members
            analysis_results = self.run_static_analysis_FE(ponding_load)
            
            # Check for convergence            
            if not self.include_ponding_effect:
                break
            
            sum_of_force = 0
            sum_of_diff  = 0
            for i in range(self.number_of_joist_spaces+1):
                for j in range(3*self.num_ele_secondary+1):            
                    sum_of_force += abs(ponding_load[i,j])
                    sum_of_diff  += abs(ponding_load_last[i,j]-ponding_load[i,j])

            print('Iteration %02i, Total Fluid Load: %.5f' % (iteration,sum_of_force))

            if sum_of_force == 0:
                if sum_of_diff <= self.tol and iteration > 0:
                    # Converged
                    break
            else:
                if sum_of_diff/sum_of_force <= self.tol:
                    # Converged
                    break

            if iteration == self.MAX_ITER-1:
                print('The maximum number iterations has been reached without convergence')
                return None
        
        if self.run_factored_analysis_after_ponding_analysis:
            analysis_results = self.run_static_analysis_FE(ponding_load,use_additional_load_factors=True)
            
        return analysis_results
        
    def run_static_analysis_FE(self,ponding_load,use_additional_load_factors=False):
    
        # Initilize Results
        results = AnalysisResults()
        results.secondary_member_deflection = np.zeros((self.number_of_joist_spaces+1,3*self.num_ele_secondary+1))
        results.primary_member_deflection_B = np.zeros((self.number_of_joist_spaces+1,1))
        results.primary_member_deflection_C = np.zeros((self.number_of_joist_spaces+1,1))
        results.primary_member_ponding_load_B = np.zeros((self.number_of_joist_spaces+1,1))
        results.primary_member_ponding_load_C = np.zeros((self.number_of_joist_spaces+1,1))
        results.secondary_member_results_AB = dict()    
        results.secondary_member_results_BC = dict()    
        results.secondary_member_results_CD = dict()    

        results.secondary_members_moment_AB = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        results.secondary_members_moment_BC = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        results.secondary_members_moment_CD = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        results.secondary_members_shear_AB = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        results.secondary_members_shear_BC = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        results.secondary_members_shear_CD = np.zeros((self.number_of_joist_spaces+1,2*self.num_ele_secondary))
        results.secondary_members_top_reaction_AB = np.zeros(self.number_of_joist_spaces+1)
        results.secondary_members_top_reaction_BC = np.zeros(self.number_of_joist_spaces+1)
        results.secondary_members_top_reaction_CD = np.zeros(self.number_of_joist_spaces+1)
        results.secondary_members_bot_reaction_AB = np.zeros(self.number_of_joist_spaces+1) 
        results.secondary_members_bot_reaction_BC = np.zeros(self.number_of_joist_spaces+1) 
        results.secondary_members_bot_reaction_CD = np.zeros(self.number_of_joist_spaces+1) 
    
        results.secondary_members_position_AB = np.zeros((2*self.num_ele_secondary,1))
        results.secondary_members_position_BC = np.zeros((2*self.num_ele_secondary,1))
        results.secondary_members_position_CD = np.zeros((2*self.num_ele_secondary,1))
        for i in range(2*self.num_ele_secondary):
            results.secondary_members_position_AB[i] = self.LAB*ceil(i/2)/self.num_ele_secondary
            results.secondary_members_position_BC[i] = self.LBC*ceil(i/2)/self.num_ele_secondary
            results.secondary_members_position_CD[i] = self.LCD*ceil(i/2)/self.num_ele_secondary
        
        results.primary_members_position = np.zeros((2*self.number_of_joist_spaces,1))
        for i in range(2*self.number_of_joist_spaces):
            results.primary_members_position[i] = self.L12*ceil(i/2)/self.number_of_joist_spaces
        
        # Secondary members between grid lines A and B
        # Build model
        secondary_member_model = FE.Model('Secondary Member AB')
        for i in range(self.num_ele_secondary+1):
            n = 'n%02i' % i
            x = (i/self.num_ele_secondary)*self.LAB
            secondary_member_model.AddNode(n,(x,0.0),('UX','UY','RZ'))
            if i == 0:
                secondary_member_model.Nodes[n].dofs['UX'].constrained = True
                secondary_member_model.Nodes[n].dofs['UY'].constrained = True
            if i == self.num_ele_secondary:
                secondary_member_model.Nodes[n].dofs['UY'].constrained = True             
            P_dead = self.dead_load_uniform*(self.L12/self.number_of_joist_spaces)*(self.LAB/self.num_ele_secondary)
            P_snow = self.snow_density*self.snow_height*(self.L12/self.number_of_joist_spaces)*(self.LAB/self.num_ele_secondary)
            if i == 0 or i == self.num_ele_secondary:
                secondary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead/2
                secondary_member_model.Nodes[n].dofs['UY'].loads['SNOW'] = -P_snow/2
            else:
                secondary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead
                secondary_member_model.Nodes[n].dofs['UY'].loads['SNOW'] = -P_snow
        for i in range(self.num_ele_secondary):
            ni = 'n%02i' % i
            nj = 'n%02i' % (i+1)
            secondary_member_model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),self.modulus_of_elasticity,self.moment_of_inertia_AB,self.cross_sectional_area)   

        # Run analyses
        for i in range(self.number_of_joist_spaces+1):
            # Apply water loads to the secondary members
            for j in range(self.num_ele_secondary + 1):
                if i == 0 or i == self.number_of_joist_spaces:
                    secondary_member_model.Nodes['n%02i' % j].dofs['UY'].loads['PONDING'] = 2*ponding_load[i,j]
                else:
                    secondary_member_model.Nodes['n%02i' % j].dofs['UY'].loads['PONDING'] = ponding_load[i,j]
                
            # Run secondary member analyses
            res = FE.LinearAnalysis(secondary_member_model)
            if use_additional_load_factors:
                res.run({
                    'DEAD':self.additional_load_factor_dead*self.alpha*self.load_factor_dead,
                    'SNOW':self.additional_load_factor_snow*self.alpha*self.load_factor_snow,
                    'PONDING':self.additional_load_factor_ponding})
            else:
                res.run({
                    'DEAD':self.alpha*self.load_factor_dead,
                    'SNOW':self.alpha*self.load_factor_snow,
                    'PONDING':1.0})

            # Get reactions for load on the primary members
            results.primary_member_ponding_load_B[i] += -secondary_member_model.Nodes['n%02i' % self.num_ele_secondary].dofs['UY'].react(res)
            
            # Get member deflections
            for j in range(self.num_ele_secondary + 1):
                results.secondary_member_deflection[i,j] = secondary_member_model.Nodes['n%02i' % j].dofs['UY'].disp(res)

            # Get member forces
            for j in range(self.num_ele_secondary):
                ele_force = secondary_member_model.Elements['e%02i'%j].force(res)/self.alpha
                results.secondary_members_moment_AB[i,2*j+0] = -ele_force.item(2)
                results.secondary_members_moment_AB[i,2*j+1] =  ele_force.item(5)
                results.secondary_members_shear_AB[i,2*j+0] =  ele_force.item(1)
                results.secondary_members_shear_AB[i,2*j+1] = -ele_force.item(4)

            # Get reactions for output
            results.secondary_members_top_reaction_AB[i] = secondary_member_model.Nodes['n00'].dofs['UY'].react(res)/self.alpha
            results.secondary_members_bot_reaction_AB[i] = secondary_member_model.Nodes['n%02i' % self.num_ele_secondary].dofs['UY'].react(res)/self.alpha          
                  
            # Save results
            results.secondary_member_results_AB[i] = res        
        
        # Secondary members between grid lines B and C
        # Build model
        secondary_member_model = FE.Model('Secondary Member BC')
        for i in range(self.num_ele_secondary+1):
            n = 'n%02i' % i
            x = (i/self.num_ele_secondary)*self.LBC
            secondary_member_model.AddNode(n,(x,0.0),('UX','UY','RZ'))
            if i == 0:
                secondary_member_model.Nodes[n].dofs['UX'].constrained = True
                secondary_member_model.Nodes[n].dofs['UY'].constrained = True
            if i == self.num_ele_secondary:
                secondary_member_model.Nodes[n].dofs['UY'].constrained = True             
            P_dead = self.dead_load_uniform*(self.L12/self.number_of_joist_spaces)*(self.LBC/self.num_ele_secondary)
            P_snow = self.snow_density*self.snow_height*(self.L12/self.number_of_joist_spaces)*(self.LBC/self.num_ele_secondary)
            if i == 0 or i == self.num_ele_secondary:
                secondary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead/2
                secondary_member_model.Nodes[n].dofs['UY'].loads['SNOW'] = -P_snow/2
            else:
                secondary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead
                secondary_member_model.Nodes[n].dofs['UY'].loads['SNOW'] = -P_snow
        for i in range(self.num_ele_secondary):
            ni = 'n%02i' % i
            nj = 'n%02i' % (i+1)
            secondary_member_model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),self.modulus_of_elasticity,self.moment_of_inertia_BC,self.cross_sectional_area)   

        # Run analyses
        for i in range(self.number_of_joist_spaces+1):
            # Apply water loads to the secondary members
            for j in range(self.num_ele_secondary + 1):
                if i == 0 or i == self.number_of_joist_spaces:
                    secondary_member_model.Nodes['n%02i' % j].dofs['UY'].loads['PONDING'] = 2*ponding_load[i,j+self.num_ele_secondary]
                else:
                    secondary_member_model.Nodes['n%02i' % j].dofs['UY'].loads['PONDING'] = ponding_load[i,j+self.num_ele_secondary]
                
            # Run secondary member analyses
            res = FE.LinearAnalysis(secondary_member_model)
            if use_additional_load_factors:
                res.run({
                    'DEAD':self.additional_load_factor_dead*self.alpha*self.load_factor_dead,
                    'SNOW':self.additional_load_factor_snow*self.alpha*self.load_factor_snow,
                    'PONDING':self.additional_load_factor_ponding})
            else:
                res.run({
                    'DEAD':self.alpha*self.load_factor_dead,
                    'SNOW':self.alpha*self.load_factor_snow,
                    'PONDING':1.0})

            # Get reactions for load on the primary members
            results.primary_member_ponding_load_B[i] += -secondary_member_model.Nodes['n00'].dofs['UY'].react(res)
            results.primary_member_ponding_load_C[i] += -secondary_member_model.Nodes['n%02i' % self.num_ele_secondary].dofs['UY'].react(res)
            
            # Get member deflections
            for j in range(self.num_ele_secondary + 1):
                results.secondary_member_deflection[i,j+self.num_ele_secondary] = secondary_member_model.Nodes['n%02i' % j].dofs['UY'].disp(res)

            # Get member forces
            for j in range(self.num_ele_secondary):
                ele_force = secondary_member_model.Elements['e%02i'%j].force(res)/self.alpha
                results.secondary_members_moment_BC[i,2*j+0] = -ele_force.item(2)
                results.secondary_members_moment_BC[i,2*j+1] =  ele_force.item(5)
                results.secondary_members_shear_BC[i,2*j+0] =  ele_force.item(1)
                results.secondary_members_shear_BC[i,2*j+1] = -ele_force.item(4)

            # Get reactions for output
            results.secondary_members_top_reaction_BC[i] = secondary_member_model.Nodes['n00'].dofs['UY'].react(res)/self.alpha
            results.secondary_members_bot_reaction_BC[i] = secondary_member_model.Nodes['n%02i' % self.num_ele_secondary].dofs['UY'].react(res)/self.alpha          
                  
            # Save results
            results.secondary_member_results_BC[i] = res  
        
        # Secondary members between grid lines C and D
        # Build model
        secondary_member_model = FE.Model('Secondary Member CD')
        for i in range(self.num_ele_secondary+1):
            n = 'n%02i' % i
            x = (i/self.num_ele_secondary)*self.LCD
            secondary_member_model.AddNode(n,(x,0.0),('UX','UY','RZ'))
            if i == 0:
                secondary_member_model.Nodes[n].dofs['UX'].constrained = True
                secondary_member_model.Nodes[n].dofs['UY'].constrained = True
            if i == self.num_ele_secondary:
                secondary_member_model.Nodes[n].dofs['UY'].constrained = True             
            P_dead = self.dead_load_uniform*(self.L12/self.number_of_joist_spaces)*(self.LCD/self.num_ele_secondary)
            P_snow = self.snow_density*self.snow_height*(self.L12/self.number_of_joist_spaces)*(self.LCD/self.num_ele_secondary)
            if i == 0 or i == self.num_ele_secondary:
                secondary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead/2
                secondary_member_model.Nodes[n].dofs['UY'].loads['SNOW'] = -P_snow/2
            else:
                secondary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead
                secondary_member_model.Nodes[n].dofs['UY'].loads['SNOW'] = -P_snow
        for i in range(self.num_ele_secondary):
            ni = 'n%02i' % i
            nj = 'n%02i' % (i+1)
            secondary_member_model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),self.modulus_of_elasticity,self.moment_of_inertia_CD,self.cross_sectional_area)   

        # Run analyses
        for i in range(self.number_of_joist_spaces+1):
            # Apply water loads to the secondary members
            for j in range(self.num_ele_secondary + 1):
                if i == 0 or i == self.number_of_joist_spaces:
                    secondary_member_model.Nodes['n%02i' % j].dofs['UY'].loads['PONDING'] = 2*ponding_load[i,j+2*self.num_ele_secondary]
                else:
                    secondary_member_model.Nodes['n%02i' % j].dofs['UY'].loads['PONDING'] = ponding_load[i,j+2*self.num_ele_secondary]
                
            # Run secondary member analyses
            res = FE.LinearAnalysis(secondary_member_model)
            if use_additional_load_factors:
                res.run({
                    'DEAD':self.additional_load_factor_dead*self.alpha*self.load_factor_dead,
                    'SNOW':self.additional_load_factor_snow*self.alpha*self.load_factor_snow,
                    'PONDING':self.additional_load_factor_ponding})
            else:
                res.run({
                    'DEAD':self.alpha*self.load_factor_dead,
                    'SNOW':self.alpha*self.load_factor_snow,
                    'PONDING':1.0})

            # Get reactions for load on the primary members
            results.primary_member_ponding_load_C[i] += -secondary_member_model.Nodes['n00'].dofs['UY'].react(res)
            
            # Get member deflections
            for j in range(self.num_ele_secondary + 1):
                results.secondary_member_deflection[i,j+2*self.num_ele_secondary] = secondary_member_model.Nodes['n%02i' % j].dofs['UY'].disp(res)

            # Get member forces
            for j in range(self.num_ele_secondary):
                ele_force = secondary_member_model.Elements['e%02i'%j].force(res)/self.alpha
                results.secondary_members_moment_CD[i,2*j+0] = -ele_force.item(2)
                results.secondary_members_moment_CD[i,2*j+1] =  ele_force.item(5)
                results.secondary_members_shear_CD[i,2*j+0] =  ele_force.item(1)
                results.secondary_members_shear_CD[i,2*j+1] = -ele_force.item(4)

            # Get reactions for output
            results.secondary_members_top_reaction_CD[i] = secondary_member_model.Nodes['n00'].dofs['UY'].react(res)/self.alpha
            results.secondary_members_bot_reaction_CD[i] = secondary_member_model.Nodes['n%02i' % self.num_ele_secondary].dofs['UY'].react(res)/self.alpha          
                  
            # Save results
            results.secondary_member_results_CD[i] = res 
    
        # Primary member on grid line B
        # Build Model
        primary_member_model = FE.Model('Primary Member B')
        for i in range(self.number_of_joist_spaces+1):
            n = 'n%02i' % i
            x = (i/self.number_of_joist_spaces)*self.L12
            primary_member_model.AddNode(n,(x,0.0),('UX','UY','RZ'))
            if i == 0:
                primary_member_model.Nodes[n].dofs['UX'].constrained = True
                primary_member_model.Nodes[n].dofs['UY'].constrained = True
            if i == self.number_of_joist_spaces:
                primary_member_model.Nodes[n].dofs['UY'].constrained = True               
            # Only apply self-weight or girder, other loads will come from the reaction of the secondary members
            P_dead = self.dead_load_on_primary_member_B*(self.L12/self.number_of_joist_spaces)
            if i == 0 or i == self.number_of_joist_spaces:
                primary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead/2
            else:
                primary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead
        for i in range(self.number_of_joist_spaces):
            ni = 'n%02i' % i
            nj = 'n%02i' % (i+1)
            primary_member_model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),self.modulus_of_elasticity,self.moment_of_inertia_B,self.cross_sectional_area)
        
        # Apply water loads to the primary members
        for i in range(self.number_of_joist_spaces+1):
            primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = results.primary_member_ponding_load_B[i]
    
        # Run analyses
        res = FE.LinearAnalysis(primary_member_model)
        if use_additional_load_factors:
            res.run({
                'DEAD':self.additional_load_factor_dead*self.alpha*self.load_factor_dead,
                'PONDING':1.0})
        else:
            res.run({
                'DEAD':self.alpha*self.load_factor_dead,
                'PONDING':1.0})
        
        # Get member defelctions
        for i in range(self.number_of_joist_spaces+1):
            results.primary_member_deflection_B[i] = primary_member_model.Nodes['n%02i' % i].dofs['UY'].disp(res)
            
        # Compute Internal Forces
        results.primary_member_moment_B = np.zeros(2*self.number_of_joist_spaces)
        results.primary_member_shear_B = np.zeros(2*self.number_of_joist_spaces)
        for i in range(self.number_of_joist_spaces):
            ele_force = primary_member_model.Elements['e%02i'%i].force(res)/self.alpha
            results.primary_member_moment_B[2*i+0] = -ele_force.item(2)
            results.primary_member_moment_B[2*i+1] =  ele_force.item(5)
            results.primary_member_shear_B[2*i+0] =  ele_force.item(1)
            results.primary_member_shear_B[2*i+1] = -ele_force.item(4)
        
        # Save results
        results.primary_member_results_B = res   
        
    
        # Primary member on grid line C
        # Build Model
        primary_member_model = FE.Model('Primary Member C')
        for i in range(self.number_of_joist_spaces+1):
            n = 'n%02i' % i
            x = (i/self.number_of_joist_spaces)*self.L12
            primary_member_model.AddNode(n,(x,0.0),('UX','UY','RZ'))
            if i == 0:
                primary_member_model.Nodes[n].dofs['UX'].constrained = True
                primary_member_model.Nodes[n].dofs['UY'].constrained = True
            if i == self.number_of_joist_spaces:
                primary_member_model.Nodes[n].dofs['UY'].constrained = True               
            # Only apply self-weight or girder, other loads will come from the reaction of the secondary members
            P_dead = self.dead_load_on_primary_member_C*(self.L12/self.number_of_joist_spaces)
            if i == 0 or i == self.number_of_joist_spaces:
                primary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead/2
            else:
                primary_member_model.Nodes[n].dofs['UY'].loads['DEAD'] = -P_dead
        for i in range(self.number_of_joist_spaces):
            ni = 'n%02i' % i
            nj = 'n%02i' % (i+1)
            primary_member_model.AddElement('e%02i'%i,'ElasticBeam2d',(ni,nj),self.modulus_of_elasticity,self.moment_of_inertia_C,self.cross_sectional_area)
        
        # Apply water loads to the primary members
        for i in range(self.number_of_joist_spaces+1):
            primary_member_model.Nodes['n%02i' % i].dofs['UY'].loads['PONDING'] = results.primary_member_ponding_load_C[i]
    
        # Run analyses
        res = FE.LinearAnalysis(primary_member_model)
        if use_additional_load_factors:
            res.run({
                'DEAD':self.additional_load_factor_dead*self.alpha*self.load_factor_dead,
                'PONDING':1.0})
        else:
            res.run({
                'DEAD':self.alpha*self.load_factor_dead,
                'PONDING':1.0})
        
        # Get member defelctions
        for i in range(self.number_of_joist_spaces+1):
            results.primary_member_deflection_C[i] = primary_member_model.Nodes['n%02i' % i].dofs['UY'].disp(res)
            
        # Compute Internal Forces
        results.primary_member_moment_C = np.zeros(2*self.number_of_joist_spaces)
        results.primary_member_shear_C = np.zeros(2*self.number_of_joist_spaces)
        for i in range(self.number_of_joist_spaces):
            ele_force = primary_member_model.Elements['e%02i'%i].force(res)/self.alpha
            results.primary_member_moment_C[2*i+0] = -ele_force.item(2)
            results.primary_member_moment_C[2*i+1] =  ele_force.item(5)
            results.primary_member_shear_C[2*i+0] =  ele_force.item(1)
            results.primary_member_shear_C[2*i+1] = -ele_force.item(4)
        
        # Save results
        results.primary_member_results_C = res  
    
                
        # Compute total roof deflection
        results.total_deflection = np.zeros((self.number_of_joist_spaces+1,3*self.num_ele_secondary+1))
        for i in range(self.number_of_joist_spaces+1):
            for j in range(self.num_ele_secondary+1):            
                x = j/self.num_ele_secondary
                # Between grid lines A and B
                primary_member_deflection = x*results.primary_member_deflection_B[i]
                results.total_deflection[i,j] = primary_member_deflection + results.secondary_member_deflection[i,j]
        
                # Between grid lines B and C
                primary_member_deflection = results.primary_member_deflection_B[i] + x*(results.primary_member_deflection_C[i] - results.primary_member_deflection_B[i])
                results.total_deflection[i,j+self.num_ele_secondary] = primary_member_deflection + results.secondary_member_deflection[i,j+self.num_ele_secondary]
                
                # Between grid lines C and D
                primary_member_deflection = (1-x)*results.primary_member_deflection_C[i]
                results.total_deflection[i,j+2*self.num_ele_secondary] = primary_member_deflection + results.secondary_member_deflection[i,j+2*self.num_ele_secondary]
        
        return results

          
def run_example():

    # Define units
    inch = 1.0
    kip = 1.0
    lb  = kip/1000.0
    ft  = 12.0*inch
    in_per_ft = inch/ft
    ksi = kip/inch**2
    plf = lb/ft
    psf = lb/ft**2
    pcf = psf/ft
    kipft = kip*ft

    # Define general parameters
    zw = 10
    strength_type = 'ASD'
    joist_span_ft = 60
    joist_girder_span_ft = 50
    joist_girder_depth_in = 48
    num_spaces = 8
    slope = 1/48

    '''
    Designation   Total Load    Defl. Limit Load    Weight
                    (plf)           (plf)           (plf)
       28K6          145              74              8.9
       28K7          162              82              9.2
       30K7          174              94              9.6
       28K8          179              89              9.8
       30K8          192             103              10
       28K9          195              97             10.5
       30K9          209             112             10.6
       28K10         232             114             11.8
       30K10         249             132             11.9
       30K11         285             150             13.3
      32LH06         294             169              14
       28K12         301             147             14.5
       30K12         324             170              15
      32LH07         329             189              16
      40LH08         348             275              16
      32LH08         357             205              17
      36LH08         361             241              18
      44LH09         443             409              19
    '''

    joist_designation_AB = '44LH09'
    joist_designation_BC = '30K11'
    joist_designation_CD = '30K11'
    joist_AB = get_joist_data(joist_designation_AB,joist_span_ft)
    joist_BC = get_joist_data(joist_designation_BC,joist_span_ft)
    joist_CD = get_joist_data(joist_designation_CD,joist_span_ft)

    # Define joists and joist girders
    joists_AB = OpenWebSteelJoist(strength_type,joist_span_ft,joist_AB.total_load('ASD'),joist_AB.deflection_limit_load())
    joists_BC = OpenWebSteelJoist(strength_type,joist_span_ft,joist_BC.total_load('ASD'),joist_BC.deflection_limit_load())
    joists_CD = OpenWebSteelJoist(strength_type,joist_span_ft,joist_CD.total_load('ASD'),joist_CD.deflection_limit_load())
    joist_girder_B = JoistGirder(strength_type,joist_girder_span_ft,joist_girder_depth_in,num_spaces,13.2)
    joist_girder_C = JoistGirder(strength_type,joist_girder_span_ft,joist_girder_depth_in,num_spaces,13.2)

    # Define class input
    input = {
        'L12': joist_girder_span_ft*ft,
        'LAB': joist_span_ft*ft,
        'LBC': joist_span_ft*ft,
        'LCD': joist_span_ft*ft,
        'number_of_joist_spaces': num_spaces,
        'dead_load_uniform': 15*psf,
        'snow_density': 15.30*pcf,
        'snow_height': 0,
        'water_density': 62.4*pcf,
        'dead_load_on_primary_member_B': 0*plf,
        'dead_load_on_primary_member_C': 0*plf,
        'alpha': 1.0,
        'load_factor_dead':    1.0,
        'load_factor_ponding': 1.0,
        'load_factor_snow':    0.0,
        'zA':  slope*0*joist_span_ft*ft,
        'zB':  slope*1*joist_span_ft*ft,
        'zC':  slope*2*joist_span_ft*ft,
        'zD':  slope*3*joist_span_ft*ft,
        'camber_of_primary_member_B': 0.0*inch,
        'camber_of_primary_member_C': 0.0*inch,
        'camber_of_secondary_members_AB': 0.0*inch,
        'camber_of_secondary_members_BC': 0.0*inch,
        'camber_of_secondary_members_CD': 0.0*inch,
        'modulus_of_elasticity': 29000*ksi,
        'cross_sectional_area': 100*inch**2,
        '_moment_of_inertia_AB': joists_AB.moment_of_inertia()/1.15,
        '_moment_of_inertia_BC': joists_BC.moment_of_inertia()/1.15,
        '_moment_of_inertia_CD': joists_CD.moment_of_inertia()/1.15,
        '_moment_of_inertia_B': joist_girder_B.moment_of_inertia()/1.15,
        '_moment_of_inertia_C': joist_girder_C.moment_of_inertia()/1.15,
        'analsis_engine': 'FE',
        'num_ele_secondary': 20,
        'include_ponding_effect': True,
    }
    print(input)

    bay = WarehouseBay(**input)
    results = bay.Run_Analysis(zw)
    
    from pprint import pprint
    #pprint(vars(bay))
    #pprint(vars(results))
    
    
    print('Utilization Ratios')
    
    utilization_ratio_joist_girder_B = 0
    utilization_ratio_joist_girder_B = max(0,
        max(joist_girder_B.moment_strength_ratio(results.primary_members_position,
                                                 results.primary_member_moment_B,
                                                 length_units='in',moment_units='kip-in')))
    utilization_ratio_joist_girder_B = max(0,
        max(joist_girder_B.shear_strength_ratio(results.primary_members_position,
                                                results.primary_member_shear_B,
                                                length_units='in',force_units='kip')))
    print(f'Joist Girder B: {utilization_ratio_joist_girder_B:.3f}')
    
    utilization_ratio_joist_girder_C = 0
    utilization_ratio_joist_girder_C = max(0,
        max(joist_girder_C.moment_strength_ratio(results.primary_members_position,
                                                 results.primary_member_moment_C,
                                                 length_units='in',moment_units='kip-in')))
    utilization_ratio_joist_girder_C = max(0,
        max(joist_girder_C.shear_strength_ratio(results.primary_members_position,
                                                results.primary_member_shear_C,
                                                length_units='in',force_units='kip')))
    print(f'Joist Girder C: {utilization_ratio_joist_girder_C:.3f}')
    
    utilization_ratio_joists_AB = 0
    for i in range(num_spaces+1):
        utilization_ratio_joists_AB = max(0,
            max(joists_AB.moment_strength_ratio(results.secondary_members_position_AB,
                                                results.secondary_members_moment_AB[i],
                                                length_units='in',moment_units='kip-in')))
        utilization_ratio_joists_AB = max(0,
            max(joists_AB.shear_strength_ratio(results.secondary_members_position_AB,
                                               results.secondary_members_shear_AB[i],
                                               length_units='in',force_units='kip')))
    print(f'Joists AB:      {utilization_ratio_joists_AB:.3f}')

    utilization_ratio_joists_BC = 0
    for i in range(num_spaces+1):
        utilization_ratio_joists_BC = max(0,
            max(joists_BC.moment_strength_ratio(results.secondary_members_position_BC,
                                                results.secondary_members_moment_BC[i],
                                                length_units='in',moment_units='kip-in')))
        utilization_ratio_joists_BC = max(0,
            max(joists_BC.shear_strength_ratio(results.secondary_members_position_BC,
                                               results.secondary_members_shear_BC[i],
                                               length_units='in',force_units='kip')))
    print(f'Joists BC:      {utilization_ratio_joists_BC:.3f}')
    
    utilization_ratio_joists_CD = 0
    for i in range(num_spaces+1):
        utilization_ratio_joists_CD = max(0,
            max(joists_CD.moment_strength_ratio(results.secondary_members_position_CD,
                                                results.secondary_members_moment_CD[i],
                                                length_units='in',moment_units='kip-in')))
        utilization_ratio_joists_CD = max(0,
            max(joists_CD.shear_strength_ratio(results.secondary_members_position_CD,
                                               results.secondary_members_shear_CD[i],
                                               length_units='in',force_units='kip')))
    print(f'Joists CD:      {utilization_ratio_joists_CD:.3f}')
    
    return


if __name__ == "__main__":
    run_example()
