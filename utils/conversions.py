def O2_PAL_to_mol(PAL):
    '''Converts O2 PAL to mol'''
    return PAL * 5.2e21 / 29
    # PAL * 5.1e14 * 1e8 / (9.8 * 32)

def mol_to_O2_PAL(mol):
    '''Converts O2 moles to PAL'''
    return mol * 29 / 5.2e21