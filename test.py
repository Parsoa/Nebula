#import rapidjson as json
import json
import time

import jellyfish

from pympler import asizeof

import plotly.offline as psklks

class Kir(object):

    class Khar(object):

        def kos(self):
            self.koon()

        def fuck(self):
            print('khar')

    def koon(self):
        a = self.Khar()
        a.kos()

    def fuck(self):
        print('kir')

if __name__ == '__main__':
    time.sleep(50)
    d = {}
    for i in range(0, 90000000):
        d[str(i)] = 'sdjkjskdjksjdkjskdjakjdlkajsdkajsldkjaslkdjalksjdklasjdkajlsdkajsld' + str(i)
    print('done')
    time.sleep(10)
    with open('./test.json', 'w') as json_file:
        json.dump(d, json_file, sort_keys = True, indent = 4)

