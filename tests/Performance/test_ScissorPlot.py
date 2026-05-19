import pytest
import numpy as np
from Objects.Performance.ScissorPlot import ScissorPlot

# ==============================
# CODE VERIFICATION TESTS
# ==============================

def test_compute_required_coefs():
    # Order of magnitude
    test_Scissorplot = ScissorPlot()
    test_Scissorplot.compute_required_coefs()
    assert test_Scissorplot.depsilon_dalpha > 0.

def test_minimum_Sh_S():
    # Order of magnitude
    test_Scissorplot_0 = ScissorPlot()
    test_Scissorplot_0.compute_required_coefs()
    assert 1 >= test_Scissorplot_0.minimum_Sh_S() > 0

    # Sensitivity
    test_Scissorplot_1 = ScissorPlot()
    test_Scissorplot_1.Vh2 = 1
    test_Scissorplot_1.compute_required_coefs()
    assert test_Scissorplot_1.minimum_Sh_S() < test_Scissorplot_0.minimum_Sh_S()







    
    