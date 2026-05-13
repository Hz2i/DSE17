import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from Objects.LinkBudget.link_budget import link_budget

# ------------------------------
# CODE VERIFICATION
# ------------------------------

# null value test
def test_link_budget_zero_distance_raises():
    budget = link_budget(distance=0, frequency=31e9, data_rate_bps=0, modulation='QPSK', target_link_margin_db=0)

    with pytest.raises(ZeroDivisionError):
        budget.compute_fsp_loss()




# print(link_budget(distance=0, frequency=31e9, data_rate_bps=0, modulation='QPSK', target_link_margin_db=0))
