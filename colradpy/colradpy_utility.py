def rate_interp_parse(rate_interp_str):
    """parse the possible interperolation rates


    :param rate_interp_str: The initial populations of levels at t=0
    :type rate_interp_str: str


    :returns: bool: log rate :interp_kind str
    """
    
    if(rate_interp_str[0:3] =='log'):
        log_rate = True
    else:
        log_rate = False

    interp_kind = 'slinear'
    if('linear' in rate_interp_str):
        interp_kind = 'linear'
    if('nearest' in rate_interp_str):
        interp_kind = 'nearest'
    if('zero' in rate_interp_str):
        interp_kind = 'zero'
    if('slinear' in rate_interp_str):
        interp_kind = 'slinear'
    if('quadratic' in rate_interp_str):
        interp_kind = 'quadratic'
    if('cubic' in rate_interp_str):
        interp_kind = 'cubic'

    return log_rate, interp_kind
