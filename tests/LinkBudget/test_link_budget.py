import pytest

from Objects.LinkBudget.link_budget import link_budget

# ==============================
# CODE VERIFICATION TESTS
# ==============================

# NULL VALUE TESTS
def test_link_budget_null_values():
    """Test link_budget calculations with zero values"""
    
    # zero distance causes ZeroDivisionError in compute_fsp_loss()
    budget = link_budget(distance=0, frequency=31e9, data_rate_bps=0, modulation='QPSK', target_link_margin_db=0)
    with pytest.raises(ZeroDivisionError):
        budget.compute_fsp_loss()
    
    # zero data_rate_bps returns 0.0 for required rx power
    result = budget.compute_required_rx_power()
    assert result == 0.0
    
    # zero distance causes ZeroDivisionError in compute_required_tx_power()
    with pytest.raises(ZeroDivisionError):
        budget.compute_required_tx_power()
    
    # zero data_rate_bps returns 0.0 for rx_sensitivity
    budget2 = link_budget(distance=1000, frequency=31e9, data_rate_bps=0, modulation='QPSK', target_link_margin_db=10)
    assert budget2.rx_sensitivity_w == 0.0


# ORDER OF MAGNITUDE TESTS
def test_link_budget_order_of_magnitude():
    """Test order-of-magnitude for key outputs"""

    distance = 400_000 # 400 km @ WikiPedia for HAPS
    budget = link_budget(distance=distance, frequency=31e9, data_rate_bps=1e6, modulation='QPSK', target_link_margin_db=10)
    
    rx_sesnsitivity = budget.rx_sensitivity_w
    assert 1e-15 < rx_sesnsitivity < 1e-13

    fsp_loss = budget.compute_fsp_loss()
    assert 1e-19 < fsp_loss < 1e-17
    
    required_rx_power = budget.compute_required_rx_power()
    assert 1e-14 < required_rx_power < 1e-12
    
    required_tx_power = budget.compute_required_tx_power()
    assert 10.0 < required_tx_power < 40.0

# ASSUMED CONSTANTS CROSS-VERIFICATION TESTS
def test_link_budget_constants():
    """Test that constants used in calculations are as expected"""
    budget = link_budget(distance=400_000, frequency=31e9, data_rate_bps=1e6, modulation='QPSK', target_link_margin_db=10)
    
    # Speed of light
    assert budget.c == 3e8
    
    # SNR requirements
    assert budget.snr_requirements_db['BPSK'] == 10.0
    assert budget.snr_requirements_db['QPSK'] == 13.0
    assert budget.snr_requirements_db['16QAM'] == 20.0
    assert budget.snr_requirements_db['64QAM'] == 26.0
    
    # Bandwidth efficiency
    assert budget.bandwidth_efficiency['BPSK'] == 1.0
    assert budget.bandwidth_efficiency['QPSK'] == 2.0
    assert budget.bandwidth_efficiency['16QAM'] == 4.0
    assert budget.bandwidth_efficiency['64QAM'] == 6.0



# UNIT TESTS
def test_link_budget_class_initialization():
    """Unit test for class initialization"""
    budget = link_budget(
        distance=12_345,
        frequency=31e9,
        data_rate_bps=2e6,
        modulation='qpsk',
        target_link_margin_db=12,
    )

    assert budget.distance == 12_345
    assert budget.frequency == 31e9
    assert budget.data_rate_bps == 2e6
    assert budget.modulation == 'QPSK'
    assert budget.target_link_margin_db == 12
    assert budget.c == 3e8


def test_link_budget_class_data_structures():
    """Unit test that class dictionaries are created correctly"""
    budget = link_budget(distance=1_000)

    assert isinstance(budget.snr_requirements_db, dict)
    assert isinstance(budget.bandwidth_efficiency, dict)
    assert set(budget.snr_requirements_db.keys()) == {'BPSK', 'QPSK', '16QAM', '64QAM'}
    assert set(budget.bandwidth_efficiency.keys()) == {'BPSK', 'QPSK', '16QAM', '64QAM'}


def test_link_budget_function_return_types():
    """Unit test for function return types on valid inputs"""
    budget = link_budget(distance=20_000, frequency=31e9, data_rate_bps=1e6, modulation='QPSK', target_link_margin_db=10)

    assert isinstance(budget.compute_fsp_loss(), float)
    assert isinstance(budget.compute_required_rx_power(), float)
    assert isinstance(budget.compute_required_tx_power(), float)


def test_link_budget_invalid_modulation_raises():
    """Unit test for expected error on unsupported modulation type"""
    with pytest.raises(KeyError):
        link_budget(distance=1_000, modulation='INVALID')