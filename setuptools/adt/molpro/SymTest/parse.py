import numpy as np


# multiple blocks of data

# file = '../ananac12.res'



# def parseResult(file):
#     ''' Parses result from the output .res files'''
#     with open(file,"r") as f:
#         dat = f.read().replace("D","E")
#     print [i.strip().split() for i in dat.strip().split("\n")]
#     # dat = [float(j) for i in dat for j in i.strip().split()]#[map(float,i.strip().split()) for i in dat]
#     dat = [[float(j) for j in i.strip().split()] for i in dat]
#     return np.array(dat)



# # def parseEnrResult

# def parseNACTResult(file):
#     data = np.genfromtxt(file)
#     # this will provide nan at places of title, so from number of nan wr can track number of blocks
#     ll=  np.where(np.isnan(data).any(axis=1))[0]
#     nBlocks = len(ll)
#     dat = np.delete(data, ll,0)
#     return np.split(dat, nBlocks)


# storage = []
# file = open(file)
# dat = file.read().replace("D","E").strip().split("\n")

# dat = [i.strip() for i in dat]
# dat = [i.split() for i in dat if i]


# getChunk = lambda seq, lim: [seq[i+1:i+lim] for i in range(0, len(seq), lim)]

# mainBlock = []
# for i in range(0, len(dat), 5):
#     datBlock = dat[i+1:i+5]
#     mainBlock.append([[ float(i) for i in j] for j in datBlock])





# def parseNactResult(file, lim):
#     with open(file, 'r') as f:
#         dat = f.read().replace("D","E").strip().split("\n")
#         dat = [i.strip() for i in dat]
#         dat = [i.split() for i in dat if i]
#     ret = [[[float(i) for i in j]
#                             for j in dat[ind+1:ind+lim]]
#                                 for ind in range(0, len(dat), lim)]
#     return np.array(ret)


# def parseEnrResult(file):
#     with open(file, "r") as f:
#         dat = f.read().replace("D","E").strip().split("\n")
#     return np.array([ float(i)  for i in dat[1:]])



# # file = '../ananac12.res'
# # print parseResult(file,5)
# # file = './energy.dat'
# data = np.loadtxt(file)

# ll = [3,2]
# l = np.cumsum(ll)[:-1]

# dat = data[:,2:]

# dats = np.split(dat, l, axis=1)
# # print dats



a = np.random.rand(2,5)
print a
ll = [0, 3,2]
l = np.cumsum(ll)
aaa = np.split(a, l, axis=1)
for i in aaa:
    print i.size