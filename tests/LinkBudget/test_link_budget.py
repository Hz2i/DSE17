import pytest
from Objects.LinkBudget import link_budget

# ------------------------------
# CODE VERIFICATION
# ------------------------------

# null value test
def test_link_budget():
    assert link_budget(0, 0, 0, 0, 0) == -float('inf')

