I/O with JSON files
'''''''''''''''''''
Krang has facilities to store all the relevant information they contain in structured JSON files. To save a JSON, the
method used, visible in the Krang reference page_, is save_json. The essential structure of such file is this:

root
|_____________"ckt_name": <string>
|_____________"elements": ______"type1.name1":
|                        |          |_________"type": <string>
|                        |          |_________"name": <string>
|                        |          |_________("units"): {....}
|                        |          |_________("depends"): {....}
|                        |          |_________("topological"): {....}
|                        |______"element2":
|                                   |....
|
|_____________"settings":______"values": {....}
|                        |_____"units": {....}
|_____________"buscoords": {....}
