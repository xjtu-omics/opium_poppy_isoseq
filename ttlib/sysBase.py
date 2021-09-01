import os

def ttSystem(cmd):
    if os.system(cmd)!=0:
        print(f'{cmd} error!')
        exit()