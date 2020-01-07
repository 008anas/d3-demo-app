def get_request_ip(request):
    """
    Returns request IP address
    """
    x_forwarded_for = request.META.get('HTTP_X_FORWARDED_FOR')
    if x_forwarded_for:
        ip = x_forwarded_for.split(',')[0]
    else:
        ip = request.META.get('REMOTE_ADDR')
    return ip

def get_session_id(request):
    if not request.session.session_key:
        request.session.save()
    return request.session.session_key