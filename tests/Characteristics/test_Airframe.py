import numpy as np
from Objects.Characteristics.Airframe import wing, fuselage, empennage, nacelles


# ==============================
# CODE VERIFICATION TESTS
# ==============================


# UNIT TESTS
def test_compute_required_coefficients_wing():
    # Sensitivity test
    test_wing_0 = wing(qc_sweep=0.0)
    test_wing_0.compute_required_coefficients()
    test_wing_1 = wing(qc_sweep=0.1)
    test_wing_1.compute_required_coefficients()
    assert test_wing_1.c_sweep > test_wing_0.c_sweep

def test_compute_CL_grad_wing():
    # Order of magnitude test
    test_wing_0 = wing(A=10**6)
    test_wing_0.compute_required_coefficients()
    test_wing_0.compute_CL_grad()
    assert 0. < test_wing_0.CL_grad < 2*np.pi

    # Sensitivity test
    test_wing_1 = wing(A=10**0)
    test_wing_1.compute_required_coefficients()
    test_wing_1.compute_CL_grad()
    assert test_wing_0.CL_grad > test_wing_1.CL_grad

    # Sensitivity test
    test_wing_2 = wing(qc_sweep=np.pi/4)
    test_wing_2.compute_required_coefficients()
    test_wing_2.compute_CL_grad()
    test_wing_3 = wing(qc_sweep=0)
    test_wing_3.compute_required_coefficients()
    test_wing_3.compute_CL_grad()
    assert test_wing_2.CL_grad < test_wing_3.CL_grad

def test_compute_CL_max():
    # Order of magnitude
    test_wing_0 = wing()
    test_wing_0.compute_required_coefficients()
    test_wing_0.compute_CL_max()
    assert 1.0 < test_wing_0.CL_max < 2.0

def test_compute_Cm_ac():# Check correct!!
    # Order of magnitude
    test_wing_0 = wing()
    test_wing_0.compute_required_coefficients()
    assert -1.0 < test_wing_0.CL_max < 1.0

# def test_compute_oswald_eff():
    # static test

def test_zero_lift_drag_wing():
    # Order of magnitude test
    test_wing_0 = wing(qc_sweep=0.0)
    test_wing_0.compute_required_coefficients()
    test_wing_0.zero_lift_drag(0.07, 25, 0.1)
    assert 0. < test_wing_0.CD0 < 0.5

    # Sensitivity test
    test_wing_1 = wing(qc_sweep=0.1)
    test_wing_1.compute_required_coefficients()
    test_wing_1.zero_lift_drag(0.07, 25, 0.1)
    assert test_wing_1.CD0 < test_wing_0.CD0

def test_zero_lift_drag_fuselage():
    # Sensitivity test
    test_fuselage_0 = fuselage(D=0.6)
    test_fuselage_0.zero_lift_drag(0.07, 25)
    test_fuselage_1 = fuselage(D=0.3)
    test_fuselage_1.zero_lift_drag(0.07, 25)
    assert test_fuselage_0.CD0 > test_fuselage_1.CD0

def test_compute_required_coefficients_empennage():
    # Sensitivity test
    test_empennage_0 = empennage(qcsweep_h=0.0)
    test_empennage_0.compute_required_coefficients()
    test_empennage_1 = empennage(qcsweep_h=0.1)
    test_empennage_1.compute_required_coefficients()
    assert test_empennage_1.qc_sweep_h > test_empennage_0.qc_sweep_h

def test_compute_Cl_grad_empennage():
    # Sensitivity test
    test_wing_0 = wing(qc_sweep=0.0)
    test_wing_0.compute_required_coefficients()
    test_wing_0.compute_CL_grad()
    test_wing_1 = wing(qc_sweep=0.1)
    test_wing_1.compute_required_coefficients()
    test_wing_1.compute_CL_grad()
    assert test_wing_1.CL_grad < test_wing_0.CL_grad

def test_zero_lift_drag_empennage():
    # Sensitivity test
    test_empennage_0 = empennage(qcsweep_h=0.0)
    test_empennage_0.compute_required_coefficients()
    test_empennage_0.zero_lift_drag(0.07, 25, 0.1)
    test_empennage_1 = empennage(qcsweep_h=0.1)
    test_empennage_1.compute_required_coefficients()
    test_empennage_1.zero_lift_drag(0.07, 25, 0.1)
    assert test_empennage_1.CD0 < test_empennage_0.CD0

def test_zero_lift_drag_nacelle():
    # Order of magnitude
    test_nacelle = nacelles()
    test_nacelle.zero_lift_drag(0.07, 25)
    assert 0 < test_nacelle.CD0 < 0.5
