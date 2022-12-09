#!/usr/bin/env python3

import sys
import math
import os.path
import struct

print(' Magicavoxel object merger script v 0.6\n----------------------------')

####################
# Input parameters #
####################

# set 'if' check to True to set the parameters in code (for dev)
if False:
    inputFile = 'M:/_projects/magicaVoxel_objectMerger/Alexia.vox'
    outputFile = 'M:/_projects/magicaVoxel_objectMerger/Alexia.merged.vox'
else:
    if len(sys.argv) < 2:
        sys.exit('\t> This script requires the address for which VOX file to merge.\n' +\
                 '\t> Also, enter an address for the VOX file to output to, unless you are ok ' +\
                    'with the default output file.\n' +
                 '\t> The default output file has the same address as the input file with ' +\
                    '".merged" added.')
    inputFile = sys.argv[1]
    if len(sys.argv) > 2:
        outputFile = sys.argv[2]
    else:
        outputFile = inputFile[:len(inputFile)-4] + '.merged' + inputFile[len(inputFile)-4:]
    
####################
# Helper functions #
####################

# Dev/Dbg functions
def _showSceneGraph(id = 0, indent = ''):
    print(indent, graphNodes[id])
    for i in range(len(graphNodes[id]['childIds'])):
        _showSceneGraph(graphNodes[id]['childIds'][i], indent + '\t')
def _checkValidMatrix(m):
    for i in range(len(m)):
        if m[i] != -1 and m[i] != 0 and m[i] != 1: return False
    if abs(m[0]) + abs(m[1]) + abs(m[2]) != 1: return False
    if abs(m[3]) + abs(m[4]) + abs(m[5]) != 1: return False
    if abs(m[6]) + abs(m[7]) + abs(m[8]) != 1: return False
    if abs(m[0] == 1) and (abs(m[3]) == 1 or abs(m[6]) == 1): return False
    if abs(m[1] == 1) and (abs(m[4]) == 1 or abs(m[7]) == 1): return False
    if abs(m[2] == 1) and (abs(m[5]) == 1 or abs(m[8]) == 1): return False
    if abs(m[3] == 1) and (abs(m[0]) == 1 or abs(m[6]) == 1): return False
    if abs(m[4] == 1) and (abs(m[1]) == 1 or abs(m[7]) == 1): return False
    if abs(m[5] == 1) and (abs(m[2]) == 1 or abs(m[8]) == 1): return False
    return True
def _reviewUniqueMatrices():
    # Find all unique rotation matrices
    matrices = []
    matricesSet = set()
    matrixToByte = {}
    for i in range(256):
        m = rotationMatrixFromDataByte(i)
        if not _checkValidMatrix(m): continue
        if str(m) in matricesSet: continue
        matrices.append(m)
        matricesSet.add(str(m))
        matrixToByte[str(m)] = i
    matrices = list(matrices)
    # Find all unique rotations from the rotation matrices
    eulers = {}
    eulersSet = set()
    for i in range(len(matrices)):
        euler = _rotationMatrixToEuler(matrices[i])
        if str(euler) in eulersSet: continue
        eulers[str(matrices[i])] = euler
        eulersSet.add(str(euler))
    # Print all unique rotations
    i = 1    
    for k, v in eulers.items():
        print(i, k, v, matrixToByte[k])
        i += 1
def _rotationMatrixToEuler(r):
    sy = math.sqrt(r[0] * r[0] +  r[3] * r[3])
    singular = sy < 1e-6
    if not singular :
        x = math.atan2(r[7], r[8])
        y = math.atan2(-r[6], sy)
        z = math.atan2(r[3], r[0])
    else :
        x = math.atan2(-r[5], r[4])
        y = math.atan2(-r[6], sy)
        z = 0
    return [int(math.degrees(x)),
            int(math.degrees(y)),
            int(math.degrees(z))]

# Parse a chunk from the data
def readChunk(data, offset):
    id = struct.unpack_from('4s', data, offset)[0]
    size1 = struct.unpack_from('i', data, offset + 4)[0]
    size2 = struct.unpack_from('i', data, offset + 8)[0]
    content = data[offset + 12 : offset + 12 + size1]
    endOffset = offset + 12 + size1 + size2
    childContent = data[offset + 12 + size1 : endOffset]
    fullChunkData = data[offset : endOffset]
    return id, content, childContent, endOffset, fullChunkData

# Parse a dictionary from the data
def readDictionary(data, startOffset):
    endOffset = startOffset
    result = {}
    dictionarySize = struct.unpack_from('i', content, endOffset)[0]
    endOffset += 4
    for i in range(dictionarySize):
        keySize = struct.unpack_from('i', content, endOffset)[0]
        endOffset += 4
        key = struct.unpack_from(str(keySize) + 's', content, endOffset)[0]
        endOffset += keySize
        valueSize = struct.unpack_from('i', content, endOffset)[0]
        endOffset += 4
        value = struct.unpack_from(str(valueSize) + 's', content, endOffset)[0]
        endOffset += valueSize
        result[key] = value
    return result, endOffset, data[startOffset: endOffset]

# Follow VOX format's rules for creating a rotation matrix based on a byte of data
def rotationMatrixFromDataByte(src):
    r = [0,0,0,0,0,0,0,0,0]
    v01 = src & 3
    v23 = (src >> 2) & 3
    v4 = (src >> 4) & 1
    v5 = (src >> 5) & 1
    v6 = (src >> 6) & 1
    # Just brute forcing this
    if v01 == 0 and v4 == 0: r[0] = 1
    if v01 == 0 and v4 == 1: r[0] = -1
    if v01 == 1 and v4 == 0: r[1] = 1
    if v01 == 1 and v4 == 1: r[1] = -1
    if v01 == 2 and v4 == 0: r[2] = 1
    if v01 == 2 and v4 == 1: r[2] = -1
    if v23 == 0 and v5 == 0: r[3] = 1
    if v23 == 0 and v5 == 1: r[3] = -1
    if v23 == 1 and v5 == 0: r[4] = 1
    if v23 == 1 and v5 == 1: r[4] = -1
    if v23 == 2 and v5 == 0: r[5] = 1
    if v23 == 2 and v5 == 1: r[5] = -1
    # Cross product for the 3rd matrix row
    row3 = crossProduct(r[0 : 3], r[3 : 6])
    r[6] = abs(row3[0])
    r[7] = abs(row3[1])
    r[8] = abs(row3[2])
    # 3rd matrix row is negative if v6 is 1
    if v6 == 1:
        r[6] *= -1
        r[7] *= -1
        r[8] *= -1
    return r

# Calculate a cross product
def crossProduct(v1, v2):
    return [ v1[1] * v2[2] - v1[2] * v2[1],
             v1[2] * v2[0] - v1[0] * v2[2],
             v1[0] * v2[1] - v1[1] * v2[0] ]

# Calculate a matrix multiplication
def multiplyMatrices(m1, m2):
    return [ m1[0] * m2[0] + m1[1] * m2[3] + m1[2] * m2[6],
             m1[0] * m2[1] + m1[1] * m2[4] + m1[2] * m2[7],
             m1[0] * m2[2] + m1[1] * m2[5] + m1[2] * m2[8],
             m1[3] * m2[0] + m1[4] * m2[3] + m1[5] * m2[6],
             m1[3] * m2[1] + m1[4] * m2[4] + m1[5] * m2[7],
             m1[3] * m2[2] + m1[4] * m2[5] + m1[5] * m2[8],
             m1[6] * m2[0] + m1[7] * m2[3] + m1[8] * m2[6],
             m1[6] * m2[1] + m1[7] * m2[4] + m1[8] * m2[7],
             m1[6] * m2[2] + m1[7] * m2[5] + m1[8] * m2[8] ]

# Transform a vector by a matrix
def transformVectorByMatrix(m, v):
    return [ v[0] * m[0] + v[1] * m[1] + v[2] * m[2],
             v[0] * m[3] + v[1] * m[4] + v[2] * m[5],
             v[0] * m[6] + v[1] * m[7] + v[2] * m[8] ] + list(v[3:])

# Add two vectors
def addToVector(v1, v2):
    return [ v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2] ] + list(v2[3:])

# Fix a MagicaVoxel bug where the translation value sometimes gets saved with a small shift from
# how it's represented in the software.  The shift is based on the rotation value (which is a byte).
# MagicaVoxel manages to fix this bug by UNshifting when the file is loaded back in.
# Unfortunately, I don't know the logic of how the to determine what to unshifting by, so I have to
# brute force the unshift based on the raw rotation byte.  The values to do this were gathered with
# painstaking experimentation.
def fixMagicalVoxelBug(rotationByte, translation):
    rotationAdjustments = {
        1:   [+0, +0, +0],  2: [+0, +0, +0],  4: [+0, +0, +0],  6: [+0, +0, +0],  8: [+0, +0, +0],
        9:   [+0, +0, +0], 17: [-1, +0, +0], 18: [-1, +0, +0], 20: [-1, +0, +0], 22: [-1, +0, +0],
        24:  [-1, +0, +0], 25: [-1, +0, +0], 33: [+0, -1, +0], 34: [+0, -1, +0], 36: [+0, -1, +0],
        38:  [+0, -1, +0], 40: [+0, -1, +0], 41: [+0, -1, +0], 49: [-1, -1, +0], 50: [-1, -1, +0],        52:  [-1, -1, +0], 54: [-1, -1, +0], 56: [-1, -1, +0], 57: [-1, -1, +0], 65: [+0, +0, -1],        66:  [+0, +0, -1], 68: [+0, +0, -1], 70: [+0, +0, -1], 72: [+0, +0, -1], 73: [+0, +0, -1],        81:  [-1, +0, -1], 82: [-1, +0, -1], 84: [-1, +0, -1], 86: [-1, +0, -1], 88: [-1, +0, -1],        89:  [-1, +0, -1], 97: [+0, -1, -1], 98: [+0, -1, -1],100: [+0, -1, -1],102: [+0, -1, -1],        104: [+0, -1, -1],105: [+0, -1, -1],113: [-1, -1, -1],114: [-1, -1, -1],116: [-1, -1, -1],        118: [-1, -1, -1],120: [-1, -1, -1],121: [-1, -1, -1]
    }
    for i in rotationAdjustments:
        if rotationByte == i:
            return addToVector(rotationAdjustments[i], translation)
    return translation

########################
# Read in the VOX data #
########################

print('\t> Reading MagicaVoxel VOX file "%s".' % inputFile)

# Get the data
if not os.path.exists(inputFile):
    sys.exit('\t>\tInput file "%s" not found.' % inputFile)
file = open(inputFile, 'rb')
data = file.read()
print('\t>\tInput file loaded (%i bytes).' % len(data))
file.close()

# Declare the info to pull from the VOX data
layerVisibilities = [ True ] * 16
modelSizes = [];
models = [];
graphNodes = [];
fullChunks = [];
# We keep the raw data of these dictionaries to simplify re-creation of the root nodes in the output
rootTrnDictionaryData = []
rootTrnFrameDictionaryData = []
rootGrpDictionaryData = []

# Header parsing
header, version = struct.unpack_from('4si', data, 0)
if header != b'VOX ':
    raise Exception('Input file is an invalid magicavoxel VOX file.')
if version != 200:
    print('\t>\t!! WARNING !! Input file has VOX version %s.  This script is intended for ' + \
          'version 200 and may not work properly.' % version)

# Parse the main chunk
id, content, chunkContent, endOffset, fullChunkData = readChunk(data, 8)
if id != b'MAIN':
    raise Exception('Input file is missing the MAIN chunk.')
if len(content) > 0:
    raise Exception('Input file has an invalid MAIN chunk.  Non-zero content.')

# Parse all chunks under the Main chunk
endOffset = 0
while endOffset < len(chunkContent):
    id, content, childContent, endOffset, fullChunkData = readChunk(chunkContent, endOffset)

    # Look for layers that are set to 'not visible', then keep the chunk data to add to the output
    if id == b'LAYR':
        index = struct.unpack_from('i', content, 0)[0]
        d, offset, dData = readDictionary(content, 4)
        if b'_hidden' in d and d[b'_hidden'] == b'1':
            layerVisibilities[index] = False
        fullChunks.append(fullChunkData)

    # parse chunks that relate to this merging
    elif id == b'SIZE':
        modelSizes.append(struct.unpack_from('iii', content, 0))
    elif id == b'XYZI':
        model = [];
        voxelCount = struct.unpack_from('i', content, 0)[0]
        offset = 4
        for i in range(voxelCount):
            model.append(struct.unpack_from('BBBB', content, offset))
            offset += 4
        models.append(model)
    elif id == b'nTRN':
        result = { 'id': -1, 'nodeType': 'TRN', 'parent': -1, 'childIds': [], 'visibility': True,
                   'layerId': -1, 'frames': [] }
        result['id'] = struct.unpack_from('i', content, 0)[0]
        d, offset, dData = readDictionary(content, 4)
        if result['id'] == 0:
            rootTrnDictionaryData = dData;
        if b'_hidden' in d and d[b'_hidden'] == b'1':
            result['visibility'] = false
        result['childIds'].append(struct.unpack_from('i', content, offset)[0])
        offset += 4
        # Skip the reserved id
        offset += 4
        result['layerId'] = struct.unpack_from('i', content, offset)[0]
        offset += 4
        frameCount = struct.unpack_from('i', content, offset)[0]
        offset += 4
        for i in range(frameCount):
            d, offset, dData = readDictionary(content, offset)
            # Frame id
            if not b'_f' in d:
                d['frameId'] = -1
            else:
                d['frameId'] = int(d[b'_f'])
                del d[b'_f']
            # Translation
            if not b'_t' in d:
                d['translation'] = [0,0,0]
            else:
                d['translation'] = [int(x) for x in d[b'_t'].split(b' ')]
                del d[b'_t']
            # MagicaVoxel bug fix - See fixMagicaVoxelBug comment for details
            if b'_r' in d:
                d['translation'] = fixMagicalVoxelBug(int(d[b'_r']), d['translation'])
            # Rotation
            if not b'_r' in d:
                d['rotation'] = [1,0,0,0,1,0,0,0,1]
            else:
                d['rotation'] = rotationMatrixFromDataByte(int(d[b'_r']))
                del d[b'_r']
            result['frames'].append(d)
            if result['id'] == 0:
                rootTrnFrameDictionaryData.append(dData);
        # Make sure the frames are sorted right
        result['frames'].sort(key=lambda o : o['frameId'])
        # Add to the collection of graphNodes
        graphNodes.append(result)
    elif id == b'nGRP':
        result = { 'id': -1, 'nodeType': 'GRP', 'parent': -1, 'childIds': [], 'visibility': True }
        result['id'] = struct.unpack_from('i', content, 0)[0]
        d, offset, dData = readDictionary(content, 4)
        if result['id'] == 1:
            rootGrpDictionaryData = dData;
        childCount = struct.unpack_from('i', content, offset)[0]
        offset += 4
        for i in range(childCount):
            result['childIds'].append(struct.unpack_from('i', content, offset)[0])
            offset += 4
        # Add to the collection of graphNodes
        graphNodes.append(result)
    elif id == b'nSHP':
        result = { 'id': -1, 'nodeType': 'SHP', 'parent': -1, 'childIds': [], 'visibility': True,
                   'framesToModelIds': {} }
        result['id'] = struct.unpack_from('i', content, 0)[0]
        d, offset, dData = readDictionary(content, 4)
        childCount = struct.unpack_from('i', content, offset)[0]
        offset += 4
        for i in range(childCount):
            modelId = struct.unpack_from('i', content, offset)[0]
            offset += 4
            d, offset, dData = readDictionary(content, offset)
            frameId = -1
            if b'_f' in d:
                frameId = int(d[b'_f'])
            result['framesToModelIds'][frameId] = modelId;
        # Add to the collection of graphNodes
        graphNodes.append(result)

    # For any chunks not related to mergeObjects: just keep the chunk data to add to the output
    else:
        fullChunks.append(fullChunkData)

######################
# Parse the VOX data #
######################

print('\t> Parsing MagicaVoxel VOX data.')

# Calculate the node parents (from the given children)
for i in range(len(graphNodes)):
    for k in range(len(graphNodes[i]['childIds'])):
        childId = graphNodes[i]['childIds'][k]
        graphNodes[childId]['parent'] = i

# hide nodes with hidden layers or hidden parent nodes
def hideChildNodes(nodeId):
    childIds = graphNodes[nodeId]['childIds']
    for i in range(len(childIds)):
        graphNodes[childIds[i]]['visibility'] = False
        hideChildNodes(childIds[i])
for i in range(len(graphNodes)):
    # Hide nodes connected to hidden layers
    if graphNodes[i]['nodeType'] == 'TRN':
        layerId = graphNodes[i]['layerId']
        if layerId >= 0 and layerVisibilities[layerId] == False:
            graphNodes[i]['visibility'] = False
    # Hide child nodes of hidden nodes
    if graphNodes[i]['visibility'] == False:
        hideChildNodes(i)

# Pass layer ids from TRN nodes to their children
for i in range(len(graphNodes)):
    if graphNodes[i]['nodeType'] != 'TRN':
        continue
    layerId = graphNodes[i]['layerId']
    childIds = graphNodes[i]['childIds']
    for k in range(len(childIds)):
        graphNodes[childIds[k]]['layerId'] = layerId

# Collect object info
voxObjects = []
for i in range(len(graphNodes)):
    if graphNodes[i]['nodeType'] != 'SHP': continue
    if graphNodes[i]['layerId'] < 0: continue
    if graphNodes[i]['visibility'] == False: continue
    newVoxObject = { 'layerId': graphNodes[i]['layerId'],
                     'frames': list(graphNodes[i]['framesToModelIds'].keys()),
                      'framesToModelIds': graphNodes[i]['framesToModelIds'], 'transforms': [] }
    # Sort frames
    newVoxObject['frames'].sort()
    # Gather the transforms
    traverseNode = graphNodes[i]
    while traverseNode['parent'] != -1:
        traverseNode = graphNodes[traverseNode['parent']]
        if traverseNode['nodeType'] == 'TRN' and len(traverseNode['frames']) > 0:
            newVoxObject['transforms'].append(traverseNode['frames'])
    # Add voxobject to the list
    voxObjects.append(newVoxObject)

# Sort voxobjects by layer
voxObjects.sort(key=lambda o : o['layerId'])

# Collect keyframes
keyframes = set()
for i in range(len(voxObjects)):
    keyframes.update(voxObjects[i]['frames'])
for i in range(len(graphNodes)):
    if graphNodes[i]['nodeType'] != 'TRN': continue
    if graphNodes[i]['visibility'] == False: continue
    for k in range(len(graphNodes[i]['frames'])):
        keyframes.add(graphNodes[i]['frames'][k]['frameId'])

keyframes = list(keyframes)
keyframes.sort()

# If there are multiple keyframes, don't use the non-key frame (for objects which aren't keyframed)
if keyframes[0] == -1 and len(keyframes) > 1:
    keyframes.pop(0)

##############################################
# Merge visible models into one final object #
##############################################

print('\t> Merging MagicaVoxel models.')

def coordToString(coord):
    return str(coord[:3])

def coordToNeighborStrings(coord):
    result = []
    result.append(coordToString([ coord[0]+1, coord[1]+0, coord[2]+0 ]))
    result.append(coordToString([ coord[0]-1, coord[1]+0, coord[2]+0 ]))
    result.append(coordToString([ coord[0]+0, coord[1]+1, coord[2]+0 ]))
    result.append(coordToString([ coord[0]+0, coord[1]-1, coord[2]+0 ]))
    result.append(coordToString([ coord[0]+0, coord[1]+0, coord[2]+1 ]))
    result.append(coordToString([ coord[0]+0, coord[1]+0, coord[2]-1 ]))
    return result

newModels = []
lowest = [ 0, 0, 0 ] # used to avoid negative space
for i in range(len(keyframes)):
    newModel = []
    filledPositions = set() # Prevent adding multiple voxels to the same coordinates
    for k in range(len(voxObjects)):
        # Find which frame to use in the vox object
        frame = -1
        for m in range(len(voxObjects[k]['frames']) -1, -1, -1):
            if voxObjects[k]['frames'][m] <= keyframes[i]:
                frame = voxObjects[k]['frames'][m]
                break
        if frame == -1:
            frame = voxObjects[k]['frames'][0]
        # Get the model
        modelId = voxObjects[k]['framesToModelIds'][frame]
        model = models[modelId]
        # Center the model by its size
        sizeAdjust = modelSizes[modelId]
        sizeAdjust = [ int(-sizeAdjust[0]/2),
                       int(-sizeAdjust[1]/2),
                       int(-sizeAdjust[2]/2) ]
        addSizeAdjustToVector = lambda v : addToVector(sizeAdjust, v)
        model = list(map(addSizeAdjustToVector, model))
        # transform the model by parent transforms
        for m in range(len(voxObjects[k]['transforms'])):
            # Find which frame to use in the transform
            transformFrame = None
            for o in range(len(voxObjects[k]['transforms'][m]) -1, -1, -1):
                if voxObjects[k]['transforms'][m][o]['frameId'] <= keyframes[i]:
                    transformFrame = voxObjects[k]['transforms'][m][o]
                    break
            if transformFrame == None:
                transformFrame = voxObjects[k]['transforms'][m][0]
            # Run the transform on the model 
            # Rotate
            rotateVector = \
                lambda v : transformVectorByMatrix(transformFrame['rotation'], v)
            model = list(map(rotateVector, model))
            # Translate
            addTranslationToVector = \
                lambda v : addToVector(transformFrame['translation'], v)
            model = list(map(addTranslationToVector, model))
        # Don't add voxels at coordinates that already have voxels
        isVectorFilled = lambda v : not coordToString(v) in filledPositions
        model = list(filter(isVectorFilled, model))
        # Keep track of where voxels are added to avoid re-adding to same position
        for m in range(len(model)):
            filledPositions.add(coordToString(model[m]))
        # Add the model to the merge
        newModel += model
    # Perform a rough hull on voxels (only removes voxels that are surrounded.  Not true hulling,
    # but better than nothing)
    hullVector = lambda v : not all(a in filledPositions for a in coordToNeighborStrings(v))
    newModel = list(filter(hullVector, newModel))
    # Calculate lowest position to avoid negative space
    for k in range(len(newModel)):
        if newModel[k][0] < lowest[0]:
            lowest[0] = newModel[k][0]
        if newModel[k][1] < lowest[1]:
            lowest[1] = newModel[k][1]
        if newModel[k][2] < lowest[2]:
            lowest[2] = newModel[k][2]
    newModels.append(newModel)

# Shift all voxels to avoid negative space (and cull any that are STILL outside boundaries)
for i in range(len(newModels)):
    shiftAllByLowest = lambda v : [ v[0] - lowest[0], v[1] - lowest[1], v[2] - lowest[2], v[3] ]
    newModels[i] = list(map(shiftAllByLowest, newModels[i]))
    # get rid of any out-of-range voxels (out-of-range voxels are exceptional)
    isByteVector = \
        lambda v : v[0]>=0 and v[0]<=255 and v[1]>=0 and v[1]<=255 and v[2]>=0 and v[2]<=255
    newModels[i] = list(filter(isByteVector, newModels[i]))

# Find model size (for all frames collectively to maintain positional consistency across frames)
newModelSize = [ -9999, -9999, -9999 ]
for i in range(len(newModels)):
    for k in range(len(newModels[i])):
        if abs(newModels[i][k][0]) > newModelSize[0]:
            newModelSize[0] = abs(newModels[i][k][0])
        if abs(newModels[i][k][1]) > newModelSize[1]:
            newModelSize[1] = abs(newModels[i][k][1])
        if abs(newModels[i][k][2]) > newModelSize[2]:
            newModelSize[2] = abs(newModels[i][k][2])
    newModelSize[0] += 1
    newModelSize[1] += 1
    newModelSize[2] += 1

###################################################
# Recreate the file content with the merged model #
###################################################

print('\t> Creating output file contents.')

# Create the model chunks
for i in range(len(newModels)):
    # Create the SIZE chunk
    print('\t>\tCreating new SIZE #%i.' % i)
    data = bytearray()
    data += struct.pack('iii', newModelSize[0], newModelSize[1], newModelSize[2]);
    contentSize = len(data)
    data = struct.pack('4s i i', b'SIZE', contentSize, 0) + data
    fullChunks.append(data)

    # Create the XYZI chunk
    print('\t>\tCreating new XYZI #%i.' % i)
    data = bytearray()
    data += struct.pack('i', len(newModels[i]))
    for k in range(len(newModels[i])):
        data += struct.pack('BBBB', newModels[i][k][0], newModels[i][k][1], \
                                    newModels[i][k][2], newModels[i][k][3])
    contentSize = len(data)
    data = struct.pack('4s i i', b'XYZI', contentSize, 0) + data
    fullChunks.append(data)

# Recreate the root Scene graph nTRN chunk
print('\t>\tRecreating root nTRN.')
data = bytearray()
data += struct.pack('i', 0)
data += rootTrnDictionaryData
data += struct.pack('iiii', 1, -1, -1, len(rootTrnFrameDictionaryData))
for i in range(len(rootTrnFrameDictionaryData)):
    data += rootTrnFrameDictionaryData[i]
contentSize = len(data)
data = struct.pack('4s i i', b'nTRN', contentSize, 0) + data
fullChunks.append(data)

# Recreate the root Scene graph nGRP chunk
print('\t>\tRecreating root nGRP.')
data = bytearray()
data += struct.pack('i', 1)
data += rootGrpDictionaryData
data += struct.pack('ii', 1, 2)
contentSize = len(data)
data = struct.pack('4s i i', b'nGRP', contentSize, 0) + data
fullChunks.append(data)

# Create the new nTRN chunk
print('\t>\tCreating new nTRN.')
data = bytearray()
data += struct.pack('iiiiii', 2, 0, 3, -1, 0, 1)
data += struct.pack('ii2s', 1, 2, b'_t')
translationText = b'0 0 ' + bytes(str(int(newModelSize[2] / 2)), 'utf-8')
data += struct.pack('i' + str(len(translationText)) + 's', len(translationText), translationText)
contentSize = len(data)
data = struct.pack('4s i i', b'nTRN', contentSize, 0) + data
fullChunks.append(data)

# Create the new nSHP chunk
print('\t>\tCreating new nSHP.')
data = bytearray()
data += struct.pack('iii', 3, 0, len(newModels))
for i in range(len(newModels)):
    keyframe = bytes(str(keyframes[i]), 'utf-8')
    data += struct.pack('iii2s', i, 1, 2, b'_f')
    data += struct.pack('i' + str(len(keyframe)) + 's', len(keyframe), keyframe)
contentSize = len(data)
data = struct.pack('4s i i', b'nSHP', contentSize, 0) + data
fullChunks.append(data)

# Create the file data with all the chunks
print('\t>\tCreating VOX file content.')
data = bytearray()
for i in range(len(fullChunks)):
    data += fullChunks[i]
contentSize = len(data)
data = struct.pack('4s i 4s i i', b'VOX ', 200, b'MAIN', 0, contentSize) + data

#########################
# Write the output file #
#########################

print('\t> Writing MagicaVoxel VOX file "' + outputFile + '".')

binary_file = open(outputFile, 'wb')
binary_file.write(data)
binary_file.close()
