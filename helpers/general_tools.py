import subprocess

def executeBinary(binpath):
    out = subprocess.run('{} -h'.format(binpath), shell=True, text=True, 
            capture_output=True)
    returncode = out.returncode
    return out, returncode
