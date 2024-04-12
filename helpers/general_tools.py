import subprocess

def executeBinary(binpath):
    args = ['-h']
    first = None
    for i in range(len(args)):
        out = subprocess.run('{} {}'.format(binpath, args[i]),
                shell=True, text=True, 
                capture_output=True)
        if out.returncode == 0:
            return out, out.returncode
        if i == 0:
            first = out
    try:
        return first, first.returncode
    except ValueError:
        return None, 1
