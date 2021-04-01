import pickle

def save_object(object, filename):
    ''' Saves a python object configuration to a pickle file

    parameters
    ----------
    object (class instance): object to be saved
    filename (string): path of the file to save the object
    '''
    with open(filename, 'wb') as f:
        pickle.dump(object, f)


def load_object(path):
    '''loads python object from a pickle file

    Parameters
    ----------
    path (string): path of the file to load the object
    '''
    with open(path, 'rb') as f:
        object = pickle.load(f)
    return object