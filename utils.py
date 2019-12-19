from astropy.io import ascii

def get_boundaries(field, id_xmin, id_xmax, id_ymin, id_ymax):
    '''
    Gets boundaries for GALFIT based on field and ID #'s of objects
    Parameters:
    -----------
    field : string
        GOODS_South or GOODS_North
    id_*** : float
        ID number of objects at boundaries of fitting region
        
    '''
    
    
    if field == 'GOODS_South':
        catalog = '/Users/rosaliaobrien/research/website/catalogs/goodss_3dhst.v4.4.cat'
        
    if field == 'GOODS_North':
        catalog = '/Users/rosaliaobrien/research/website/catalogs/goodsn_3dhst.v4.4.cat'
        
    cat = ascii.read(catalog)
    x_min = cat['x']['id'==id_xmin].shape[1]-100
    x_max = cat['x']['id'==id_xmax].shape[1]+100
    y_min = cat['y']['id'==id_ymin].shape[1]-100
    y_max = cat['y']['id'==id_ymax].shape[1]+100
    print(x_min, x_max, y_min, y_max)