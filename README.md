# MagicaVoxel - Object merger
A python 3 script to create a modified version of a VOX file.  The modification is to merge visible objects in all layers down to a single object while maintaining animation key-frames.

---

Example 1:

`python objectMerger.py Alexia.vox`

Example 2:

`python objectMerger.py Alexia.vox Alexia.merged.vox`

---

The second parameter (`Alexia.merged.vox` in example 2) is optional.  It represents what file to write the output to.  The second parameter in example 2 is the same as the default, which example 1 uses implicitly.

`Note:` Dragging the "Alecia.vox" file onto the script has the same effect as example 1.

---

These examples will load the MagicaVoxel VOX file "Alexia.vox" (which is included in this repository for demonstration purposes and has many objects in many layers and a few animation key-frames).

Once "Alexia.vox" is loaded, a new file "Alexia.merged.vox" will be created.  "Alexia.merged.vox" will be a duplicate of "Alexia.vox", except that all objects in all visible layers will be flattened into a single object on the first layer.  All animation key-frames will be maintained during the merge.

