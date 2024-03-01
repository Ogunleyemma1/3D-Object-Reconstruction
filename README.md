# 3D-Object-Reconstruction
Projective and direct Euclidean Reconstruction from 2D Images
With knowledge of the relative orientation, spatial object coordinates can be triangulated from corresponding image points. If the parameters of the interior orientation are unknown, then only a projective reconstruction is possible. Using at least five control points, this intermediate result can be transformed quite simply into a Euclidean reconstruction.
Projective reconstruction:
Since the manual matching of image points is quite laborious and boring, a text file bh.dat with many homologous image points is made available for the image pair showing the bust of BEETHOVEN.
a) Read the homologous image coordinates 12↔xxin the format (x1, y1, x2, y2), e.g. with
fh = fopen('bh.dat', 'r'); A = fscanf(fh, '%f%f%f%f', [4 inf]); fclose(fh);
x1 = A(1:2, :); x2 = A(3:4, :);
and use your function from exercise 4 in order to determine the relative orientation of the images with the fundamental matrix F.
b) Implement a new function, which defines two corresponding projection matrices PN and P’ by means of F.
c) Realize a function for the linear triangulation of projective object points XP1 and try to visualize the computed spatial object coordinates, e.g. using
figure; scatter3(X(1,:), X(2,:), X(3,:), 10, 'filled'); axis square; view(32, 75);
Direct Euclidean reconstruction:
a) Read the control point information from the provided file pp.dat in the format (x1, y1, x2, y2, XE, YE, ZE) and triangulate projective object points XP2 from the five homologous image points 12↔xx using the already computed projection matrices PN and P’.
b) Extend your algorithm from exercise 2 for the planar 2D homography to a spatial 3D homography H. Determine the spatial transformation of the five projective object points XP2 to the corresponding Euclidean object points XE.
c) Apply this transformation H to all object points of your projective reconstruction XP1 and visualize the result of the Euclidean reconstruction spatially.
