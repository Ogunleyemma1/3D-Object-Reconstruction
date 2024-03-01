function ProjectiveReconstruction

%________________PART 1____________________________________

%Get the Image Cordinate Points
[X1, X2] = GetImageCordinates();

%Obtain the Fundamental Matrix
fMatrix = GetFundamentalMatrix(X1,X2);

%Perform Singularity Constraint on the Fundamental Matrix
FundamentalMatrix = SingularityConstraint(fMatrix);

%Obtain Projection Matrices for the two Image Points
[P1, P2] = CreateProjectionMatrix(FundamentalMatrix);

%Obtain the Object Coordinates by Linear Triangulation
ObjectPoints1 = LinearTriangulation(X1, X2, P1, P2);

%Plot Object Points
PlotObjectPoint(ObjectPoints1);




%________________PART 2____________________________



%Reading Coordinate points from Data provided
[C1, C2, C3] = ReadControlPoint();

%Obtaining Object Points by Linear Triangulation
ObjectPoint2 = LinearTriangulation(C1, C2, P1, P2);
disp(ObjectPoint2);

%Performing 3D Homography on the Object Points
H = Get3DHomography(ObjectPoint2, C3);

% Applying the homography to Object Points 1
EuclieanObjectPoint = GetEuclieanObjectPoint(H, ObjectPoints1);

%Visualize ObjectPoint 1
PlotObjectPoint(EuclieanObjectPoint);

end



%______________PART 1 METHODS__________________________________________

function [ImageCoordinate1, ImageCoordinate2] = GetImageCordinates()

fh = fopen("bh.dat","r");
A = fscanf(fh, '%f%f%f%f', [4 inf]);
fclose(fh);
ImageCoordinate1 = A(1:2, :);
ImageCoordinate2 = A(3:4, :);
ImageCoordinate1(3,:) = 1;
ImageCoordinate2(3,:) = 1;
end

function fMatrix = GetFundamentalMatrix(ImageCoordinate1, ImageCoordinate2)

%Performing translation, scaling and conditioning

%Translation
Image1Midpoint = [mean(ImageCoordinate1(1,:)); mean(ImageCoordinate1(2,:)); mean(ImageCoordinate1(3,:))];
Image2Midpoint = [mean(ImageCoordinate2(1,:)); mean(ImageCoordinate2(2,:)); mean(ImageCoordinate2(3,:))];

Image1Translated = [(ImageCoordinate1(1,:)-Image1Midpoint(1, 1)); (ImageCoordinate1(2,:)-Image1Midpoint(2, 1)); (ImageCoordinate1(3,:)-Image1Midpoint(3, 1))];
Image2Translated = [(ImageCoordinate2(1,:)-Image2Midpoint(1, 1)); (ImageCoordinate2(2,:)-Image2Midpoint(2, 1)); (ImageCoordinate2(3,:)-Image2Midpoint(3, 1))];

%scaling
Image1Scaled = [mean(abs(Image1Translated(1,:))); mean(abs(Image1Translated(2,:))); mean(abs(Image1Translated(3,:)))];
Image2Scaled = [mean(abs(Image2Translated(1,:))); mean(abs(Image2Translated(2,:))); mean(abs(Image2Translated(3,:)))];

%Translation Matrix
Image1TranslatedMatrix = [1, 0, -Image1Midpoint(1,1); 0, 1, -Image1Midpoint(2,1); 0, 0, 1];
Image2TranslatedMatrix = [1, 0, -Image2Midpoint(1,1); 0, 1, -Image2Midpoint(2,1); 0, 0, 1];

%Scaled Matrix
Image1ScaledMatrix = [1/Image1Scaled(1,1), 0, 0; 0, 1/Image1Scaled(2,1), 0; 0, 0, 1];
Image2ScaledMatrix = [1/Image2Scaled(1,1), 0, 0; 0, 1/Image2Scaled(2,1), 0; 0, 0, 1];

%Image 1 and 2 Transformation Matrix
Image1TransformationMatrix = Image1ScaledMatrix * Image1TranslatedMatrix;
Image2TransformationMatrix = Image2ScaledMatrix * Image2TranslatedMatrix;

%Image 1 and 2 Conditioned Coordinates
Image1Conditioned = [Image1Translated(1,:)*Image1TransformationMatrix(1,1); Image1Translated(2,:)*Image1TransformationMatrix(2,2); Image1Translated(3,:)];
Image2Conditioned = [Image2Translated(1,:)*Image2TransformationMatrix(1,1); Image2Translated(2,:)*Image2TransformationMatrix(2,2); Image2Translated(3,:)];

%Formulating the Design Matrix

A = zeros(size(Image1Conditioned, 2), 9);

for i = 1:size(Image1Conditioned,2)
    
    A(i, :) = [Image1Conditioned(1, i)*Image2Conditioned(1, i), Image1Conditioned(2, i)*Image2Conditioned(1, i), Image2Conditioned(1, i), Image1Conditioned(1, i)*Image2Conditioned(2, i), Image1Conditioned(2, i)*Image2Conditioned(2, i), Image2Conditioned(2, i), Image1Conditioned(1, i), Image1Conditioned(2, i), 1];
    
end

%Applying SVD
[U, D, V] = svd(A);

H = reshape(V(:,end), 3, 3)';

fMatrix = Image2TransformationMatrix' * H * Image1TransformationMatrix; %Computing the reverse conditioning

fMatrix = fMatrix(:,:)/fMatrix(end,end);  %Normalizing
end

function FundamentalMatrix = SingularityConstraint(fMatrix)
%Enforcing Singularity Constraint
if det(fMatrix) == 0

    FundamentalMatrix = fmatrix;
    disp("The Fundamental matrix is:");
    disp(FundamentalMatrix);

else

    [U, D, V] = svd(fMatrix);
    
    % replacing the last diagonal element of D with 0

    D(end, end) = 0;
    FundamentalMatrix = U * D * V';
    disp("The Fundamental matrix is:");
    disp(FundamentalMatrix);
end 
end

function [PN, P] = CreateProjectionMatrix(FundamentalMatrix)
%Projection Matrix Perpendicular to the Object for Image 1
PN = [eye(3) zeros(3,1)]; 

%Calulating the epipoles of Image 2
[U, D, V] = svd(FundamentalMatrix);

Image2Ep = U(:,end);
disp("The Epipole of Image 2 is:");
disp(Image2Ep);

%Obtaining the skew-symmetric Matrix
aX = [0, -Image2Ep(3,1), Image2Ep(2,1);
    Image2Ep(3,1), 0, -Image2Ep(1,1);
    -Image2Ep(2,1), Image2Ep(1,1), 0];

%Obtaining the projection Matrix for Image 2
P = (aX * FundamentalMatrix) + [Image2Ep, Image2Ep, Image2Ep];
P = [P, Image2Ep];

disp("Image1 Projective Matrix is:");
disp(PN);
disp("Image2 Projection Matrix is:");
disp(P);
end


function ObjectPoints = LinearTriangulation(ImageCoordinate1, ImageCoordinate2, PN, P)

% Creating the Object Point Design Matrix
A = [];

% Initialize ObjectPoints to an empty matrix
ObjectPoints = [];

for i = 1:size(ImageCoordinate1, 2)
    A = [ImageCoordinate1(1,i)*PN(3,:)-PN(1,:)
         ImageCoordinate1(2,i)*PN(3,:)-PN(2,:)
         ImageCoordinate2(1,i)*P(3,:)-P(1,:)
         ImageCoordinate2(2,i)*P(3,:)-P(2,:)];

    % Performing Singular Value Decomposition
    [Uo, Do, Vo] = svd(A);

    % Extracting the 3D coordinates from the last column of Vo
    % Accumulate the results in ObjectPoints
    ObjectPoints = [ObjectPoints, Vo(:, size(Vo,2)) ./ Vo(4, size(Vo,2))];
    end
end


function PlotObjectPoint(ObjectPoints)

figure; scatter3(ObjectPoints(1,:), ObjectPoints(2,:), ObjectPoints(3,:), 10, "filled");
axis square; view(32, 75);

end



%________________PART 2 METHODS__________
function [X1, X2, X3] = ReadControlPoint()
fh = fopen("pp.dat","r");
A = fscanf(fh, '%f%f%f%f', [7 inf]);
fclose(fh);
X1 = A(1:2, :);
X2 = A(3:4, :);
X3 = A(5:7, :);
X1(3,:) = 1;
X2(3,:) = 1;
X3(4, :) = 1;

end

function Homography = Get3DHomography(ImageCoordinate1, ImageCoordinate2)

%Performing translation, scaling and conditioning

%Translation
Image1Midpoint = [mean(ImageCoordinate1(1,:)); mean(ImageCoordinate1(2,:)); mean(ImageCoordinate1(3,:))];
Image2Midpoint = [mean(ImageCoordinate2(1,:)); mean(ImageCoordinate2(2,:)); mean(ImageCoordinate2(3,:))];

Image1Translated = [(ImageCoordinate1(1,:)-Image1Midpoint(1, 1)); (ImageCoordinate1(2,:)-Image1Midpoint(2, 1)); (ImageCoordinate1(3,:)-Image1Midpoint(3, 1))];
Image2Translated = [(ImageCoordinate2(1,:)-Image2Midpoint(1, 1)); (ImageCoordinate2(2,:)-Image2Midpoint(2, 1)); (ImageCoordinate2(3,:)-Image2Midpoint(3, 1))];

%scaling
Image1Scaled = [mean(abs(Image1Translated(1,:))); mean(abs(Image1Translated(2,:))); mean(abs(Image1Translated(3,:)))];
Image2Scaled = [mean(abs(Image2Translated(1,:))); mean(abs(Image2Translated(2,:))); mean(abs(Image2Translated(3,:)))];

%Translation Matrix
Image1TranslatedMatrix = [1, 0, 0, -Image1Midpoint(1,1); 0, 1, 0, -Image1Midpoint(2,1); 0, 0, 1, -Image1Midpoint(3,1); 0, 0, 0, 1];
Image2TranslatedMatrix = [1, 0, 0, -Image2Midpoint(1,1); 0, 1, 0, -Image2Midpoint(2,1); 0, 0, 1, -Image2Midpoint(3,1); 0, 0, 0, 1];

%Scaled Matrix
Image1ScaledMatrix = [1/Image1Scaled(1,1), 0, 0, 0; 0, 1/Image1Scaled(2,1), 0, 0; 0, 0, 1/Image1Scaled(3,1), 0; 0, 0, 0, 1];
Image2ScaledMatrix = [1/Image2Scaled(1,1), 0, 0, 0; 0, 1/Image2Scaled(2,1), 0, 0; 0, 0, 1/Image2Scaled(3,1), 0; 0, 0, 0, 1];

%Image 1 and 2 Transformation Matrix
Image1TransformationMatrix = Image1ScaledMatrix * Image1TranslatedMatrix;
Image2TransformationMatrix = Image2ScaledMatrix * Image2TranslatedMatrix;
disp(Image1TransformationMatrix);
disp(Image2TransformationMatrix);

%Image 1 and 2 Conditioned Coordinates
Image1Conditioned = [Image1Translated(1,:)*Image1TransformationMatrix(1,1); Image1Translated(2,:)*Image1TransformationMatrix(2,2); Image1Translated(3,:)*Image1TransformationMatrix(3,3)];
Image2Conditioned = [Image2Translated(1,:)*Image2TransformationMatrix(1,1); Image2Translated(2,:)*Image2TransformationMatrix(2,2); Image2Translated(3,:)*Image2TransformationMatrix(3,3)];

%Formulating the Design Matrix

A = [];

for i = 1:size(Image1Conditioned,2)
    
    A = [A; -Image1Conditioned(1, i), -Image1Conditioned(2, i), -Image1Conditioned(3, i), -1, zeros(1, 4), zeros(1,4), Image1Conditioned(1,i)*Image2Conditioned(1, i), Image1Conditioned(2,i)*Image2Conditioned(1, i), Image1Conditioned(3,i)*Image2Conditioned(1, i), Image2Conditioned(1, i);
        zeros(1,4), -Image1Conditioned(1, i), -Image1Conditioned(2, i), -Image1Conditioned(3, i), -1, zeros(1,4), Image1Conditioned(1,i)*Image2Conditioned(2, i), Image1Conditioned(2,i)*Image2Conditioned(2, i), Image1Conditioned(3,i)*Image2Conditioned(2, i), Image2Conditioned(2, i);
        zeros(1,4), zeros(1,4), -Image1Conditioned(1, i), -Image1Conditioned(2, i), -Image1Conditioned(3, i), -1, Image1Conditioned(1,i)*Image2Conditioned(3, i), Image1Conditioned(2,i)*Image2Conditioned(3, i), Image1Conditioned(3,i)*Image2Conditioned(3, i), Image2Conditioned(3, i)];
end
disp(A);


%Applying SVD
[U, D, V] = svd(A);
disp(V);

H = reshape(V(:,end), 4, 4)';


Homography = inv(Image2TransformationMatrix) * H * Image1TransformationMatrix; %Computing the reverse conditioning

Homography = Homography(:,:)/Homography(end,end);  %Normalizing
disp("The 3D Homography Matrix is:");
disp(Homography);
end


function EuclideanObjectPoints = GetEuclieanObjectPoint(Homography, ObjectPoint)

EuclideanObjectPoints = [];

Points = Homography*ObjectPoint;

for i = 1:size(ObjectPoint, 2)

    EuclideanObjectPoints = [EuclideanObjectPoints, Points(:, i)./Points(4,i)];
end
end


