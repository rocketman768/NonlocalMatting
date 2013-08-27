# Nonlocal Matting

This is the research code used by me, Philip G. Lee, in my nonlocal matting
research. Please see the paper here:
[Nonlocal matting](http://users.eecs.northwestern.edu/~pgl622/files/NonlocalMatting_Lee_2011.pdf)

## License

All of the code is licensed under the GPLv3. The images in gt/ and images/ are from
[www.alphamatting.com](http://www.alphamatting.com) under fair use for
research. The images in scribs/ are licensed under the WTFPL-2. For more
detail, see the **COPYRIGHT** file.

## Purpose

The purpose of this code was to do research on incorporating Nonlocal methods
into alpha matting. The results of the original code are published in
<i>P. Lee, Y Wu.
["Nonlocal matting"](http://users.eecs.northwestern.edu/~pgl622/files/NonlocalMatting_Lee_2011.pdf).
Computer Vision and Pattern Recognition 2011.</i>

## Notes

This code is meant for **RESEARCH ONLY**. It is not production quality or
optimized or anything like that. It is written to be easily modified.

### Normalized Images

Most of this code assumes that we are working with images whose color channels
are normalized to the range [0,1]. For example:
    I = double(imread('image.png'))./255
is the normalization for an 8-bit image.

### Memory Limitations

The code directly solves sparse matrices in memory, limiting the image
resolution to something like 128x128, since it will exhaust GB of RAM.
