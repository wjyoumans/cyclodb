import os, sys
import json
from sage.all import *


def format(line):
    js = json.loads(line)

    keys = ["conductor", "degree", "discriminant_bits", "signature", "norm_relation", "precision"]
    for k in keys:
        v = js[k]
        js[k] = {"plain": v, "filter": v, "display_raw": f"{v}", "display": f"\\\\({v}\\\\)"}

    # no processing but use floats for filter
    #keys = ["discriminant", "h", "h_minus", "h_plus"]
    #for k in keys:
    #    v = js[k]
    #    js[k] = {"plain": v, "filter": float(v), "display": f"\\({v}\\)"}
   

    # real numbers (input as strings)
    keys = ["regulator", "residue"]
    for k in keys:
        v = js[k]

        temp = str(RR(v)).split('e')
        if len(temp) == 1:
            s = f"{temp[0]}"
        else:
            s = f"{temp[0]} \\\\times 10^{{{temp[1].strip('+')}}}"
        js[k] = {"plain": v, "filter": float(RR(v)), "display_raw": s, "display": f"\\\\({s}\\\\)"}

    # large ints
    keys = ["discriminant", "h", "h_minus", "h_plus"]
    for k in keys:
        v = js[k]

        temp = str(RR(v)).split('e')
        if len(temp) == 1:
            s = f"{v}"
        else:
            s = f"{temp[0]} \\\\times 10^{{{temp[1].strip('+')}}}"
        js[k] = {"plain": v, "filter": int(v), "display_raw": s, "display": f"\\\\({s}\\\\)"}

    # polynomials
    keys = ["polynomial"]
    Z = PolynomialRing(ZZ, "x")
    for k in keys:
        v = js[k]
        s = f"{latex(Z(v))}"
        js[k] = {"plain": v, "filter": v, "display_raw": s, "display": f"\\\\({s}\\\\)"}
    
    # finite abelian groups

    keys = ["class_group", "galois_group"]
    orders = ["h", "degree"]
    for k, ko in zip(keys, orders):
        v = js[k]

        s = ""
        if len(v) == 0:
            s = "C_1"
        else:
            set_v = list(reversed(sorted(list(set(v)))))
            for i in range(len(set_v)-1):
                m = set_v[i]
                exp = v.count(m)
                if exp == 1:
                    s += f"C_{{{m}}} \\\\times "
                else:
                    s += f"C_{{{m}}}^{{{exp}}} \\\\times "

            m = set_v[-1]
            exp = v.count(m)
            if exp == 1:
                s += f"C_{{{m}}}"
            else:
                s += f"C_{{{m}}}^{{{exp}}}"

        js[k] = {"plain": v, "filter": js[ko]["filter"], "display_raw": s, "display": f"\\\\({s}\\\\)"}
                
    return json.dumps(js)

def main():
    #DIR = os.path.dirname(os.path.abspath(__file__))

    json_out = "cyclodata.json"
    DIR = sys.argv[1]

    total = 0
    err = 0
    good = 0

    temp = []
    for f in os.listdir(DIR):
        #if total > 10:
        #    break
        total += 1

        if ".swp" in f:
            print('swp')
            continue
        
        with open(os.path.join(DIR, f), "r") as fp:
            lines = fp.readlines()
            if not "sys" in lines[-1]:
                err += 1
                continue
           
            for line in lines:
                if 'stack' in line:
                    err +=1
                    break
                if line[0] == '{':
                    good += 1
                    temp.append(format(line))
                    break


    print(f"Errors: {err}")
    print(f"Good: {good}")
    print(f"Total: {total}")
    assert(total == err + good)

    s = "{ \"data\": [" + ",".join(temp).strip('\n') + "]}"
    js = json.loads(s)

    with open(json_out, "w") as fp:
        json.dump(js, fp)

if __name__ == "__main__":
    main()
