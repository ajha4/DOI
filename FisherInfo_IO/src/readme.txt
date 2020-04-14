Each exp file is named as Expx_y.bin
x is the no of voxels ( x = 20 for low scatter exp ( mus = 1 cm-1) and = 40 for medium scatter (mus = 4 cm-1))
y is the experiment no
y = 11 : homogeneous media ( mua = 0.01 cm -1)
y = 12 : homogeneous media ( mua = 0.1 cm -1)
y = 2 : three different tissue types scattered randomly. 2nd tiss type has mua = 1 cm-1 and 3rd has mus = 0.1 cm-1;
y = 3 : spherical blob with high absorption ( mua = 1 cm-1) at the center
y = 4 : spherical blob with very low scatter and absorption ( mua = 0.001 cm-1) at the center

FIM files:
All files have absorption co. of 0.25 mm
_depth.dat: Depth varies from 1.5 mm to 17.5 mm (in res of 1 mm). Sc. co. varies from 0.25 cm-1
to 2.5 cm-1 (Spacing of 0.25 cm). FIM arranged as depth varies first across column, and sc. co will
vary across row. Radius of signal = 1.5 mm
_size.dat: Depth varies from 1.5 mm to 17.5 mm (spacing of 1 mm). Size
(radius) varies from 0.5 mm to 2.5 mm(Spacing of 0.5 mm). FIM arranged as size varies first
across column, and depth will vary across row. Sc. co of  signal = 1 cm-1
_offcenter.dat: Only depth varies from 1.5 mm to 17.5 mm (in res of 1 mm). Sc.
co. is 1cm-1, and signal radius = 1.5 mm
