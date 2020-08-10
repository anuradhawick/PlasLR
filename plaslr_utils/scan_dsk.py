import h5py
from multiprocessing import Pool
import sys
import os
import argparse
import logging

logger = logging.getLogger('PlasLR')

def write(args):
    key, dsk_output, output = args
    
    cmd = f"""h5dump -d /dsk/solid/"""+ str(key) +f""" -y --ddl "{dsk_output}" | awk '{{if($1=="{{" || $1 == "}}," || $1=="}}") {{ printf "" }} else if (match($1, /,$/) || $1 == "") {{ printf $1}}  else print $1}}' > "{output}/{key}.chunk" """        
    os.system(cmd)
    
    logger.debug("Finished chunk - " + str(key))


def scan_dsk(dsk_output, threads, output):
    logger.info("Reading DSK result")
    dsk_h5 = h5py.File(dsk_output, "r")
    keys = list(map(str, dsk_h5["dsk/solid"].keys()))
    dsk_h5.close()

    p = Pool(threads)

    p.map(write, [(key, dsk_output, output) for key in keys])
    p.close()
    p.join()

    logger.info("Gathering DSK result")
    
    cmd = """cat "{0}"/*.chunk > "{0}/15mersCounts" """.format(output)
    os.system(cmd)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Filter reads and ground truth with threshold 1000bp""")

    parser.add_argument('--dsk', '-i',
                        help="DSK input",
                        type=str,
                        required=True)
    parser.add_argument('--output', '-o',
                        help="Output input",
                        type=str,
                        required=True)   
    parser.add_argument('--threads', '-t',
                        help="Thread count",
                        type=int,
                        required=True)                        
    

    args = parser.parse_args()
    dsk_output = args.dsk
    output = args.output
    threads = args.threads

    scan_dsk(dsk_output, threads, output)

    
