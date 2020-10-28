%Project 2 400073796
clear all
close all
clc

% number of rows of your reference image. Use a small number if your computer is slow
row = 100;
img_ori = im2double(imresize(rgb2gray(imread('Capture2.PNG')), [row NaN]));
figure
imshow(img_ori)
title('Reference Image')


% the linear motion of a camera by len pixels horizontally
len = round(0.1 * size(img_ori,1));

for i = 1:size(img_ori,1)
    for j = 1:size(img_ori,2)-len+1
        img_motion(i, j) = mean(img_ori(i, j:j+len-1));
    end
end
figure
imshow(img_motion)
title('Motion Blur Image')


motion_matrix = zeros(numel(img_motion), numel(img_ori));
order = reshape((1:numel(img_motion))', size(img_motion));
for i = 1:size(img_ori,1)
    fprintf('building motion matrix...%.1f%%\n', i/size(img_ori,1)*100);
    for j = 1:size(img_ori,2)-len+1
        img_temp = zeros(size(img_ori));
        img_temp(i, j:j+len-1) = 1.0/len;
        motion_matrix(order(i, j), :) = reshape(img_temp, 1, []);
    end
end


img_matrix_blur = reshape(motion_matrix * reshape(img_ori, [], 1), size(img_motion));
if max(max(abs(img_matrix_blur - img_motion))) > 1e-10
    error('wrong motion matrix');
end


boundary_matrix = zeros((len-1)*size(img_ori,1), numel(img_ori));
boundary_counter = 1;
for i = 1:size(img_ori,1)
    fprintf('building boundary matrix...%.1f%%\n', i/size(img_ori,1)*100);
    for j = 1:len-1
        img_temp = zeros(size(img_ori));
        img_temp(i, j) = 1.0;
        boundary_matrix(boundary_counter, :) = reshape(img_temp, 1, []);
        boundary_vector(boundary_counter, 1) = img_ori(i, j);
        boundary_counter = boundary_counter + 1;
    end
end

%Out of focus blur matrix
%A= zeros(100);
% for i=1:100
%     for j=1:100
%         if i==j
%             A(i,j)=1/2;
%         
%         end
%         if i-1==j
%             A(i,j)=1/8;
%         end
%         if i+1==j
%             A(i,j)=1/8;
%         end
%     end
% end


% linear equations: Ax = b
A = [motion_matrix; boundary_matrix];
b = [reshape(img_motion, [], 1); boundary_vector];
img_deblur = reshape(Solving_Linear_Equations_with_LU_decomposition(A, b), size(img_ori));
  
  
if max(max(abs(img_deblur - img_ori))) > 1e-10
    error('wrong deblur image');
end
figure
imshow(img_deblur)
title('Deblur Image')
disp('Done!')


%Displaying Solved Image
%TEST CASES
 %t= [1 1 -1; 1 -2 3; 2 3 1];
% 
 %p=[4; -6; 7];
% 
 %Solving_Linear_Equations_with_LU_decomposition(t, p)


%LU Decomp Function
function x = Solving_Linear_Equations_with_LU_decomposition(A,B)
   
   %step 1.

    %[L,U] = lu(A); %finding lower and upper matrix *uncomment for
    %efficiency
    
    U=zeros(size(A,1)); %setting zero matrix for U and L
    L=zeros(size(A,2));
    
    for i=1:size(A,1)
        L(i,i)=1; %make diagonal line of 1's
        U(1,i)=A(1,i); %clone first row
    end
    dimension = size(A,1); 
    for i=2:dimension % row start at 2
        for j=1:dimension % col start at 2
            %Lower matrix logic
            for q=1:i-1 %range of q is one less the number of row we are on
                s1=0;
                if q==1 %first q is always 0
                    s1=0;
                else
                for p=1:q-1
                    s1=s1+L(i,p)*U(p,q);
                end
                end
                L(i,q)=(A(i,q)-s1)/U(q,q); 
            end
            %upper matrix logic
             for q=i:dimension
                 s2=0;
               for p=1:i-1
                   s2=s2+L(i,p)*U(p,q);
               end
               U(i,q)=A(i,q)-s2; %upper matrix doesn't depend on numbers of lower
               
             end
        end
    end
   
    %step 2.
    %Y=linsolve(L,B); %solving for Y in LY=B *uncomment for efficiency
    
    
    for i = 1:size(L,1) %traversing from row
        Y(i)= B(i); %intializing that Y(i) = B(i) at start
        j=i; 
        while j>1 %traversing aslong as column is shorter than size
            Y(i)=Y(i)-L(i,j-1)*Y(j-1);
            j=j-1;
        end
        Y(i)=Y(i)/L(i,i);
    end
   
    %step 3.
    %x=linsolve(U,Y); %solve for x in UX=Y *uncomment for efficiency
    
     for i = size(U,1):-1:1 %traversing from row
        X(i)= Y(i); %intializing that Y(i) = B(i) at start
        j=i; 
        while j<size(U,1) %traversing aslong as column is shorter than size
            X(i)=X(i)-U(i,j+1)*X(j+1);
            j=j+1;
        end
        X(i)=X(i)/U(i,i);
     end
     X=X'

    
   
end