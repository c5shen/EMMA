'''
Created on 6.27.2022 by Chengze Shen

for the sake of sub-modules
'''

from configs import Configs

def notifyError(location):
    print('Encountered an error at {}, check {}'.format(
        location, Configs.error_path))
    exit(1)
