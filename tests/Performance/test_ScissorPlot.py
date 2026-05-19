import pytest
import numpy as np
from Objects.Performance.ScissorPlot import ScissorPlot
from Objects.Characteristics.Airframe import wing, empennage, fuselage
from Objects.Characteristics.ReferenceGeometries import airfoil_e387, airfoil_NACA0012

# ==============================
# CODE VERIFICATION TESTS
# ==============================


def make_scissor_plot():
    return ScissorPlot(
        wing=wing(airfoil=airfoil_e387()),
        empennage=empennage(airfoilh=airfoil_NACA0012(), airfoilv=airfoil_NACA0012()),
        fuselage=fuselage(),
    )

def test_compute_required_coefs():
    # Order of magnitude
    test_Scissorplot = make_scissor_plot()
    test_Scissorplot.compute_required_coefs()
    assert test_Scissorplot.depsilon_dalpha > 0.

def test_minimum_Sh_S():
    # Order of magnitude
    test_Scissorplot_0 = make_scissor_plot()
    test_Scissorplot_0.compute_required_coefs()
    assert 1 >= test_Scissorplot_0.minimum_Sh_S() > 0

    # Sensitivity
    test_Scissorplot_1 = make_scissor_plot()
    test_Scissorplot_1.Vh2 = 1
    test_Scissorplot_1.compute_required_coefs()
    assert test_Scissorplot_1.minimum_Sh_S() < test_Scissorplot_0.minimum_Sh_S()







    
    