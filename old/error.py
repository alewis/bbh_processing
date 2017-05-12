import os
def error(msg):
    os.sys.stderr.write("#### ERROR ####\n");
    os.sys.stderr.write(msg)
    os.sys.stderr.write('\n')
    os.sys.exit(-1)
