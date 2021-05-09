import random


def read():
    r = []
    with open('pi_pos_all.fasta', 'r') as f:
        for l in f:
            if l.startswith('>'):
                continue
            r.append(l.strip())
    return r

def main(target):
    r = read()
    n = len(r)
    m = len(r[0])
    for i in range(target):
        s = ''
        while len(s) < m:
            seq = random.randint(0, n - 1)
            auglen = random.randint(3, m - len(s)) // 3 * 3
            s += r[seq][len(s):len(s) + auglen]
        print('>OEAVGG{:06d}'.format(i))
        print(s)


if __name__ == '__main__':
    target = 10000
    main(target)
