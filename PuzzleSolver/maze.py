import os

input_file = 'cqi.txt'

if not os.path.exists(input_file):
    print(f"Error: The input file '{input_file}' does not exist in the current directory.")
    print("Current directory contents:")
    print(os.listdir())
    exit(1)

with open(input_file, 'r', encoding='utf-8') as file:
    content = file.read()

lines = content.strip().split('\n')

data = []
home = []
sinks = []
connections = []

min_x = 0
max_x = 0
max_y = 0

for line in lines:
    parts = line.split()
    if len(parts) == 3:
        uchar, x, y = parts
        xc = int(x)
        yc = int(y)
        if xc > max_x: max_x=xc
        if yc > max_y: max_y=yc
        if(uchar=="*"):
            home = [uchar,xc,yc]
        data.append([uchar,xc,yc])

maze = sorted(data, key=lambda item: (-item[2], item[1]))

conn = {
    "╠": [True, True, True, False],
    "╚": [True, True, False, False],
    "║": [True, False, True, False],
    "╩": [True, True, False, True],
    "╣": [True, False, True, True],
    "╝": [True, False, False, True],
    "╦": [False, True, True, True],
    "╔": [False, True, True, False],
    "╗": [False, False, True, True],
    "═": [False, True, False, True]
}


def xyToID(x,y):
    return (y * (max_x - min_x)) + x

def idToXY(gid):
    y = gid // max_x
    if y != 0 and gid % y == 0 and gid % max_x == 0:
        x = gid // y
    else:
        x = gid % max_x
    return [x,y]

for i, item in enumerate(maze):
    if i >= 0:
        if item[1] != maze[i - 1][1] + 1:
            print("", end="\n" if item[1] == max_x else " ")
    print(f"{item[0]}", end="\n" if item[1] == max_x else "")
   
visited = [xyToID(home[1],home[2])]

def checkConn(i,ii):
    if i < 0 or ii < 0:
        return False
    else:
        if [i,ii] not in connections and [ii,i] not in connections:
            connections.append([i,ii])
            if ii not in visited:
                visited.append(ii)
                tt,ts = idToXY(ii)
            return True

def search_maze(x, y):
    for item in maze:
        if item[1] == x and item[2] == y:
            return item
    return None


   
for i, item in enumerate(visited):
    
    x, y = idToXY(visited[i])
    now = search_maze(x,y)
    nid = xyToID(x,y)
    u=xyToID(x,y+1)
    r=xyToID(x+1,y)
    d=xyToID(x,y-1)
    l=xyToID(x-1,y)
    ids = [u,r,d,l]

    up = search_maze(x,y+1)
    right = search_maze(x+1,y)
    down = search_maze(x,y-1)
    left = search_maze(x-1,y)
    neighbors = [up,right,down,left]

    sequence = [2, 3, 0, 1]
    

    def charCheck(char):
        if char=='*' or char==' ' or char.isalpha():
            return False
        else:
            return True

    def checkCheck(a,aa):
        if(a == False and aa == False):
            return False
        elif(a==aa):
            return True
        else:
            return False
   
    def checkDirection(it,other,c):
        if charCheck(it[0]):
            if charCheck(other[0]):
                g=conn[it[0]][c]
                o=other[0]
                h=conn[o][sequence[c]]
                if checkCheck(g,h):
                    if checkConn(nid,ids[c]):
                        return True
            else:
                if checkConn(nid,ids[c]):
                    return True
        else:
            if checkConn(nid,ids[c]):
                return True

    try:
        if now[0].isalpha() or now[0]=='*':
            if nid not in sinks:
                sinks.append(nid)
            checkConn(nid,u)
            checkConn(nid,r)
            checkConn(nid,d)
            checkConn(nid,l)
        elif now[0]==' ':
            pass#print("***",nid,"***")
        else:

            cc = now[0]
            cpn = 0
            for s in range(4):
                try:
                    if neighbors[s]==None:
                        pass#print("NOne")
                    else:
                        if checkDirection(now,neighbors[s],s):
                            cpn = cpn +1
                            t,tt = idToXY(ids[s])
                            if ids[s] not in visited:
                                visited.insert(i+cpn,ids[s])
                except Exception as e: pass#print(e)

    except Exception as e: pass#print(e)

print("conn",len(connections))
print("sinks:",len(sinks))

for t, itemt in enumerate(sinks):
    x, y = idToXY(sinks[t])
    now = search_maze(x,y)
    print(now," - ",sinks[t])




