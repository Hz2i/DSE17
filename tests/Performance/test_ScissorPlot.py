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
    # assert test_Scissorplot.x_ac - test_Scissorplot.wing.x_ac > 0







    
    