from mimetypes import init
import numpy as np
from Objects.Characteristics.Airframe import wing, fuselage, empennage


def test_compute_CL_grad():
    # Order of magnitude test
    test_wing_1 = wing(A=10**6)
    test_wing_1.compute_CL_grad()
    assert 0. < test_wing_1.CL_grad < 2*np.pi

    # Sensitivity test
    test_wing_2 = wing(A=10**0)
    test_wing_2.compute_CL_grad()
    assert test_wing_1.CL_grad > test_wing_2.CL_grad

    # Sensitivity test
    test_wing_3 = wing(qc_sweep=10)
    test_wing_3.compute_CL_grad()
    test_wing_4 = wing(qc_sweep=0)
    test_wing_4.compute_CL_grad()
    assert test_wing_3.CL_grad < test_wing_4.CL_grad

# def test_compute_CL_max():
    # not complete

# def test_compute_oswald_eff():
    # static test

def test_zero_lift_drag_wing():
    # Order of magnitude test
    test_wing_5 = wing(qc_sweep=10)
    test_wing_5.zero_lift_drag()
    assert 0. < test_wing_5.CD0 < 0.1

    # Sensitivity test
    test_wing_6 = wing(qc_sweep=0)
    test_wing_6.zero_lift_drag()
    assert test_wing_5.CD0 < test_wing_6.CD0

def test_zero_lift_drag_fuselage():
    # Sensitivity test
    test_fuselage_1 = fuselage(D=0.6)
    test_fuselage_1.zero_lift_drag()
    test_fuselage_2 = fuselage(D=0.3)
    test_fuselage_2.zero_lift_drag()
    assert test_fuselage_1.CD0 > test_fuselage_2.CD0

def test_zero_lift_drag_wing_empennage(): # not complete
    # Sensitivity test
    test_empennage_1 = empennage()
    test_empennage_1.zero_lift_drag()
    test_empennage_2 = empennage()
    test_empennage_2.zero_lift_drag()
    assert test_empennage_1.CD0 > test_empennage_2.CD0