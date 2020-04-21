def get_session_id(request):
    if not request.session.session_key:
        request.session.save()
    return request.session.session_key


def data_to_html_data(data, datatype, filename=None):
    from base64 import b64encode

    """Data types: zip, genbank, fasta, pdf"""
    datatype = {
        "zip": "application/zip",
        "genbank": "application/genbank",
        "fasta": "application/fasta",
        "pdf": "application/pdf",
        "xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    }.get(datatype, datatype)
    datatype = 'data:%s;' % datatype
    data64 = 'base64,%s' % b64encode(data).decode("utf-8")
    headers = ''
    if filename is not None:
        headers += "headers=filename%3D" + filename + ";"
    return datatype + headers + data64
