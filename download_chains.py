#!/usr/bin/python3
import sys
import subprocess

# scp -r -P 13900 joaoreboucas@che.cbpf.br:~/cocoa/Cocoa/projects/cs2-project/chains/ .
PATH="~/cocoa2/Cocoa/projects/cs2_project/chains/"
if __name__ == "__main__":
    match len(sys.argv):
        case 1:
            print("Downloading all chains...")
            proc = subprocess.run([
                "scp", "-r", "-P", "13900", f"joaoreboucas@che.cbpf.br:{PATH}", "./chains/"
            ])
            if proc.returncode != 0: print("ERROR downloading all chains")
        case 2:
            indices = []
            # Parse indices such as '23', '10-20', '1,2,3-20'
            for input in sys.argv[1].split(','):
                if "-" in input:
                    imin, imax = list(map(int, input.split('-')))
                    for i in range(imin, imax + 1): indices.append(i)
                else:
                    indices.append(int(input))
            print(f"Downloading chains {indices}")
            for i in indices:
                proc = subprocess.run([
                    "scp", "-r", "-P", "13900", f"joaoreboucas@che.cbpf.br:{PATH}/MCMC{i}", "./chains/"
                ])
                if proc.returncode != 0: print(f"ERROR downloading chain {i}")
        case default:
            print("ERROR: unsupported amount of arguments")
            exit(1)