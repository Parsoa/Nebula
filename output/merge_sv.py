import io
import json

tracks = {}
for i in range(0, 12):
    with open('./batch_' + str(i) + '.json') as json_file:
        batch = json.load(json_file)
        tracks.update(batch)
with open('./inversion_boundaries.json', 'w') as json_file:
    json.dump(tracks, json_file, sort_keys=True, indent=4, separators=(',', ': '))

