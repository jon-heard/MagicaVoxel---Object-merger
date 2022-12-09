# MagicaVoxel---Object-merger
Python script to merge visible objects in all layers down to a single object, retaining animation frames.

Example command-line executions:

```
python objectMerger.py Alexia.vox

python objectMerger.py Alexia.vox Alexia.merged.vox
```

Both examples above have the same effect.  The second parameter is optional, and the second parameter entered in the second example is the same as the default (used in the first example).

These examples will load the MagicaVoxel VOX file "Alexia.vox" (which is included in this repository for demonstration purposes and has many objects in many layers and a few animation key-frames).

Once "Alexia.vox" is loaded, a new file "Alexia.merged.vox" will be created.  "Alexia.merged.vox" will be a duplicate of "Alexia.vox", except that all objects in all visible layers will be flattened into a single object on the first layer.  All key-frames, both from voxels and objects transforms, will be maintained during the flattening, allowing the animation to remain intact.
