"""
This file is intended to handle the contingencies applied during the design process.
It will take the current values of the design parameters and check them against the target values.
If the current value exceeds the target value, it will raise a warning.
The contingencies have to be updated manually by the risk manager or the risk manager has to be informed of intended changes.
Ideally, the contingencies reduce to zero when the design is finalised.
Inputs: current values from design tools, target values from budgets, contingencies from risk manager
Outputs: TPMs, warning messages (subsystem, parameter, current value, target value, approach and contingency)
"""

class ContingencyManagement:
    def __init__(self, design_stage="detailed"):                                   
        """Initialise the contingency management system with predefined contingencies based on the design stage. 
        The contingencies are defined as percentages and can be adjusted as needed, 
        but that requires authorization of the risk manager.
        """
        if design_stage == "conceptual":
            self.comp_mass_contingency = 0.3 # 30% mass contingency
            self.comp_power_contingency = 0.3 # 30% power contingency
            self.comp_volume_contingency = 0.3 # 30% volume contingency
            self.comp_cost_contingency = 0.3 # 30% cost contingency
            
            self.comms_mass_contingency = 0.3 # 30% mass contingency
            self.comms_power_contingency = 0.3 # 30% power contingency
            self.comms_mass_contingency = 0.3 # 30% mass contingency
            self.comms_mass_contingency = 0.3 # 30% mass contingency
            
            self.fcs_mass_contingency = 0.3 # 30% mass contingency
            self.fcs_power_contingency = 0.3 # 30% power contingency
            self.fcs_volume_contingency = 0.3 # 30% volume contingency
            self.fcs_cost_contingency = 0.3 # 30% cost contingency
            
            self.paysys_mass_contingency = 0.3 # 30% mass contingency
            self.paysys_power_contingency = 0.3 # 30% power contingency
            self.paysys_volume_contingency = 0.3 # 30% volume contingency
            self.paysys_cost_contingency = 0.3 # 30% cost contingency
            
            self.ctrl_mass_contingency = 0.3 # 30% mass contingency
            self.ctrl_power_contingency = 0.3 # 30% power contingency
            self.ctrl_volume_contingency = 0.3 # 30% volume contingency
            self.ctrl_cost_contingency = 0.3 # 30% cost contingency
            
        elif design_stage == "preliminary":
            self.comp_mass_contingency = 0.2 # 20% mass contingency
            self.comp_power_contingency = 0.2 # 20% power contingency
            self.comp_volume_contingency = 0.2 # 20% volume contingency
            self.comp_cost_contingency = 0.2 # 20% cost contingency
            
            self.comms_mass_contingency = 0.25 # 25% mass contingency
            self.comms_power_contingency = 0.25 # 25% power contingency
            self.comms_volume_contingency = 0.25 # 25% volume contingency
            self.comms_cost_contingency = 0.25 # 25% cost contingency
            
            self.fcs_mass_contingency = 0.25 # 25% mass contingency
            self.fcs_power_contingency = 0.25 # 25% power contingency
            self.fcs_volume_contingency = 0.25 # 25% volume contingency
            self.fcs_cost_contingency = 0.25 # 25% cost contingency
            
            self.paysys_mass_contingency = 0.25 # 25% mass contingency
            self.paysys_power_contingency = 0.25 # 25% power contingency
            self.paysys_volume_contingency = 0.25 # 25% volume contingency
            self.paysys_cost_contingency = 0.25 # 25% cost contingency
            
            self.ctrl_mass_contingency = 0.2 # 20% mass contingency
            self.ctrl_power_contingency = 0.2 # 20% power contingency
            self.ctrl_volume_contingency = 0.2 # 20% volume contingency
            self.ctrl_cost_contingency = 0.2 # 20% cost contingency

        elif design_stage == "detailed":
            self.comp_mass_contingency = None #TBD 
            self.comp_power_contingency = None #TBD 
            self.comp_volume_contingency = None #TBD 
            self.comp_cost_contingency = None #TBD
            
            self.comms_mass_contingency = None #TBD
            self.comms_power_contingency = None #TBD
            self.comms_volume_contingency = None #TBD
            self.comms_cost_contingency = None #TBD 
            
            self.fcs_mass_contingency = None #TBD 
            self.fcs_power_contingency = None #TBD 
            self.fcs_volume_contingency = None #TBD
            self.fcs_cost_contingency = None #TBD
            
            self.paysys_mass_contingency = None #TBD 
            self.paysys_power_contingency = None #TBD
            self.paysys_volume_contingency = None #TBD 
            self.paysys_cost_contingency = None #TBD 
            
            self.ctrl_mass_contingency = None #TBD
            self.ctrl_power_contingency = None #TBD 
            self.ctrl_volume_contingency = None #TBD
            self.ctrl_cost_contingency = None #TBD 
        
        elif design_stage == "final":
            self.comp_mass_contingency = None #TBD
            self.comp_power_contingency = None #TBD 
            self.comp_volume_contingency = None #TBD 
            self.comp_cost_contingency = None #TBD 
            
            self.comms_mass_contingency = None #TBD 
            self.comms_power_contingency = None #TBD 
            self.comms_volume_contingency = None #TBD 
            self.comms_cost_contingency = None #TBD 

            self.fcs_mass_contingency = None #TBD 
            self.fcs_power_contingency = None #TBD 
            self.fcs_volume_contingency = None #TBD 
            self.fcs_cost_contingency = None #TBD 
            
            self.paysys_mass_contingency = None #TBD
            self.paysys_power_contingency = None #TBD
            self.paysys_volume_contingency = None #TBD 
            self.paysys_cost_contingency = None #TBD 
            
            self.ctrl_mass_contingency = None #TBD 
            self.ctrl_power_contingency = None #TBD 
            self.ctrl_volume_contingency = None #TBD 
            self.ctrl_cost_contingency = None #TBD 
        
        else:
            raise ValueError('Invalid design stage. Please choose from "conceptual", "preliminary", "detailed", or "final".')
        
        # actual values from tools
        self.comp_mass_actual = None # Has to be imported from somewhere
        self.comp_power_actual = None 
        self.comp_volume_actual = None
        self.comp_cost_actual = None

        self.comms_mass_actual = None
        self.comms_power_actual = None
        self.comms_volume_actual = None
        self.comms_cost_actual = None

        self.fcs_mass_actual = None
        self.fcs_power_actual = None
        self.fcs_volume_actual = None
        self.fcs_cost_actual = None

        self.paysys_mass_actual = None
        self.paysys_power_actual = None
        self.paysys_volume_actual = None
        self.paysys_cost_actual = None

        self.ctrl_mass_actual = None
        self.ctrl_power_actual = None
        self.ctrl_volume_actual = None
        self.ctrl_cost_actual = None

        # target values from budgets
        self.comp_mass_target = None # Has to be imported from somewhere
        self.comp_power_target = None
        self.comp_volume_target = None
        self.comp_cost_target = None

        self.comms_mass_target = None
        self.comms_power_target = None
        self.comms_volume_target = None
        self.comms_cost_target = None

        self.fcs_mass_target = None
        self.fcs_power_target = None
        self.fcs_volume_target = None
        self.fcs_cost_target = None

        self.paysys_mass_target = None
        self.paysys_power_target = None
        self.paysys_volume_target = None
        self.paysys_cost_target = None

        self.ctrl_mass_target = None
        self.ctrl_power_target = None
        self.ctrl_volume_target = None
        self.ctrl_cost_target = None

    def calculate_current_value(self, actual_value, contingency, approach="overestimate"):
        """Calculate the current value based on the actual value and the contingency percentage."""
        if approach == "underestimate":
            return actual_value / (1 + contingency)
        elif approach == "overestimate":
            return actual_value * (1 + contingency)
        else:
            raise ValueError('Invalid approach. Please choose "underestimate" or "overestimate".')

    def assess_risk(self, subsystem, parameter, current_value, target_value, approach, contingency):                        
        """
        Assess the risk by comparing the current value with the target value.
        If the current value is not compatible with the target value, a warning message is raised 
        indicating the subsystem, parameter, current value, target value, approach and contingency.
        """
        # Define general warning message format
        warning_message = (f"Warning: The conservative {current_value} is not compatible with the"
                            f" target value of {target_value} for the parameter {parameter} under the {approach} approach "
                            f"with a contingency of {contingency*100}% in the {subsystem} subsystem.")
        confirmation_message = (f"Confirmation: The conservative {current_value} is compatible with the"
                                f" target value of {target_value} for the parameter {parameter} under the {approach} approach "
                                f"with a contingency of {contingency*100}% in the {subsystem} subsystem.")
        # Compare the critical current value with the target value and raise a warning if necessary
        if approach == "underestimate" and current_value < target_value or approach == "overestimate" and current_value > target_value:
            print(warning_message)
        else:
            print(confirmation_message)
    
if __name__ == "__main__":
    # Example usage
    #contingency_manager = ContingencyManagement(design_stage="conceptual")
    #current_mass = contingency_manager.calculate_current_value(contingency_manager.comp_mass_actual, contingency_manager.comp_mass_contingency, approach="overestimate")
    #contingency_manager.assess_risk(subsystem="Computer System", parameter="Mass", current_value=current_mass, target_value=contingency_manager.comp_mass_target, approach="overestimate", contingency=contingency_manager.comp_mass_contingency)
    pass