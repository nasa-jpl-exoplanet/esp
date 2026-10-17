'''
GMR:Core code for SV views
Can only return a plt.figure()
'''

import matplotlib.pyplot as plt

def timing(dct:dict, title=''):
    '''
    GMR:Timing SV figure
    '''
    fgr = plt.figure(figsize=(12, 9))
    plt.title(title, fontsize=20)
    plt.plot(dct['time'], dct['z'], 'o')
    plt.xlabel('[MJD-UTC]', fontsize=20)
    plt.ylabel(r'Separation [R$^*$]', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    return fgr
