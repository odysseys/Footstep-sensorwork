# system test
import sys
import string
import time
import thread
import Queue
import struct
import os.path
import fnmatch
import matplotlib.pyplot as plt
import numpy as np

bufferSize = 256
# DATA_FOLDER = './copyFromPi/Footstep/'
DATA_FOLDER = './'
dataAll = []
TIME_OVERFLOW = 1099511627776
TIME_RES = 0.000015650040064103
data = ''

for (dir, dirs, files) in os.walk(DATA_FOLDER):
    for file in files:
        print(file)
        if fnmatch.fnmatch(file, '*.txt'):
            try:
                inputStr = DATA_FOLDER + file
                outStr1 = DATA_FOLDER + file + 'out.txt'
                f = open(inputStr,'rb')
                f1 = open(outStr1,'w')

                dataFile = f.read()
                data = data + dataFile
                strList = data.split("\xFF\xFF")
                strListLen = len(strList)
                if len(strList[strListLen-1]) != 5 and len(strList[strListLen-1]) != 4+5 and len(strList[strListLen-1]) != bufferSize:
                    data = strList[strListLen-1]
                    del strList[strListLen-1]
                    # print("not finished")
                else:
                    data = ''

                for vs in strList:
                    if len(vs) == 0:    
                        continue
                    if len(vs) == 5:
                        timestamp = bytearray(b"\x00\x00\x00")
                        timestamp.append(vs[0])
                        timestamp.append(vs[1])
                        timestamp.append(vs[2])
                        timestamp.append(vs[3])
                        timestamp.append(vs[4])
                        timestamp = struct.unpack("L", timestamp)[0]
                        # print(vs[0:4])
                        if timestamp < 0:
                            timestamp = math.fmod(float(timestamp), TIME_OVERFLOW) * TIME_RES
                        print(timestamp)
                        f1.write(str(timestamp) + '\n')
                        # qr.put(item = timestamp,block = False)
                        # f1.write(str(timestamp) + '\n')
                        # timestamp1 =  ord(vs[7])+ord(vs[8])*256+ord(vs[9])*256^2+ord(vs[10])*256^3+ord(vs[11])*256^4
                        # # print(vs[7:11])
                        # timestamp2 =  ord(vs[14])+ord(vs[15])*256+ord(vs[16])*256^2+ord(vs[17])*256^3+ord(vs[18])*256^4
                        # print(vs[14:18])
                    elif len(vs) == 4+5: 
                        # timestamp =  ord(vs[0])+ord(vs[1])*256+ord(vs[2])*256^2+ord(vs[3])*256^3+ord(vs[4])*256^4
                        # print(timestamp)
                        distance = vs[5:8]
                        distance = struct.unpack("f", vs[5:9])[0]
                        f1.write(str(d) + '\n')
                        # qr.put(item = distance,block = False)
                        print(distance)
                        # f1.write(str(distance) + '\n')
                    elif len(vs) == bufferSize: 
                        # print("data")
                        for i in range(0,bufferSize-1,2):
                            d = ord(vs[i])+ord(vs[i+1])*256
                            f1.write(str(d) + '\n')
                            dataAll.append(d)
                    else: 
                        print("wrong: "+ str(len(vs))) 
                        if len(vs) == bufferSize+1:
                            for i in range(1,bufferSize,2):
                                d = ord(vs[i])+ord(vs[i+1])*256
                                f1.write(str(d) + '\n')
                                dataAll.append(d)
            finally:
                f.close()
                f1.close()

            # plt.figure()
            # plt.plot(dataAll)