# MagicaVoxel - Object merger
A python 3 script to create a modified version of a VOX file.  The modification is to merge visible objects in all layers down to a single object while retaining animation frames.

Example command-line executions:

```
python objectMerger.py Alexia.vox

python objectMerger.py Alexia.vox Alexia.merged.vox
```

Both examples above have the same effect.  The second parameter is optional.  The second example's second parameter is the same as the default (which the first example uses implicitly).

These examples will load the MagicaVoxel VOX file "Alexia.vox" (which is included in this repository for demonstration purposes and has many objects in many layers and a few animation key-frames).

Once "Alexia.vox" is loaded, a new file "Alexia.merged.vox" will be created.  "Alexia.merged.vox" will be a duplicate of "Alexia.vox", except that all objects in all visible layers will be flattened into a single object on the first layer.  All key-frames, both from voxels and from objects transforms, will be maintained during the merge, allowing the animations to remain intact.
