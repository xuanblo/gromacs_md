import re 
import sys

last_position = 0
all_pos = 1
pos = 1

with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
    for line in lines:
        if "SOL" not in line:
            obj = re.match(r'\s+(\d+)[A-Z]{3}', line)
            if obj:
                if obj.group(1) == last_position:
                    pass
                else:
                    if int(last_position) + 1 == int(obj.group(1)):
                        pos += 1
                    else:
                        pos = 1
                    last_position = obj.group(1)
                    print(all_pos, pos, obj.group(1))
                    all_pos += 1