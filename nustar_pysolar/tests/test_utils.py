import nustar_pysolar.utils as utils

def test_channel_to_keV():
    """Channel to keV check."""
    assert utils.channel_to_keV(210) == 10.
    
def test_keV_to_channel():
    """keV to channel check."""
    assert utils.keV_to_channel(10) == 210
    
    
def test_convert_nustar_time():
    """Time conversion check."""
    assert utils.convert_nustar_time(80000000).isoformat() == '2012-07-14T22:13:15'
    
    