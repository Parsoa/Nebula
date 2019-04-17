import os
import sys

if __name__ == '__main__':
   
    with open(os.path.join('test1.json'), 'w') as f:
        f.write('workdir')
        print >> sys.stderr, '1'
    with open(os.path.join('/', 'test2.json'), 'w') as f:
        f.write('/\n')
        print >> sys.stderr, '2'
    with open(os.path.join('/share', 'test3.json'), 'w') as f:
        f.write('/share\n')
        print >> sys.stderr, '3'
    with open(os.path.join('/share/hormozdiarilab', 'test4.json'), 'w') as f:
        f.write('/share/hormozdiarilab\n')
        print >> sys.stderr, '4'
    with open(os.path.join('/share/hormozdiarilab/Codes', 'test5.json'), 'w') as f:
        f.write('/share/hormozdiarilab/Codes\n')
        print >> sys.stderr, '5'
    with open(os.path.join('/share/hormozdiarilab/Codes/NebulousSerendipity', 'test6.json'), 'w') as f:
        f.write('/share/hormozdiarilab/Codes/NebulousSerendipity\n')
        print >> sys.stderr, '6'
    with open(os.path.join('/share/hormozdiarilab/Codes/NebulousSerendipity/src', 'test7.json'), 'w') as f:
        f.write('/share/hormozdiarilab/Codes/NebulousSerendipity/src\n')
        print >> sys.stderr, '7'
    with open(os.path.join('/share/hormozdiarilab/Codes/NebulousSerendipity/src/python', 'test8.json'), 'w') as f:
        f.write('/share/hormozdiarilab/Codes/NebulousSerendipity/src/python\n')
        print >> sys.stderr, '8'
    with open(os.path.join(os.getcwd(), 'test9.json'), 'w') as f:
        f.write(os.getcwd() + '\n')
        print >> sys.stderr, '9'

    exit()
