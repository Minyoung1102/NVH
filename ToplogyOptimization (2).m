function xPhys = ToplogyOptimization(w,w_load,VolRatio,Filter_Image)

Active = Filter_Image;
nelx = 128;
nely = 128;
%% Parameters

indent = 2;
Ncell = 1;
rad = 0.15 *nelx;
betaMax = 10000; % Max. beta 값.
SetVol = 0.4; % Predefined field 가 차지할 volume.
volfrac = mean(Active(:))*VolRatio;
x(1:nely,1:nelx) = volfrac*ones(nely,nelx);
savefig = 0;
%% Fixed parameters

penal = 3;
rmin = 1.8;
eta = 0.5; 
Mnd2_old = 100; 
beta = 1; 
List = 0.01:0.01:1; % 적절한 threshold 를 찾기 위한 후보 리스트 작성.
E0 = 210000.; 
Emin = 1e-3;
nu = 0.3; 
FixNumber = [];
LoadBC = [];
FixNode = [];
Fix = zeros(nely,nelx);
loop = 0;
change = 1.; % Initial change 값.
% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% Filter Preparing
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);



%% Outer boundary 와 Fixed Boundary condition 설정.

for ely = 1 : nely+1
    for elx = 1 : nelx+1
        
            % Load 가 걸릴 node 좌표 정하기.
            if (ely-round((nely+2)/2))^2+(elx-round((nelx+2)/2))^2 >= (nely/2*0.92-indent)^2 && (ely-round((nely+2)/2))^2+(elx-round((nelx+2)/2))^2 <= (nely/2*0.97-indent)^2
                LoadBC(end+1,:) = [ely,elx];
            end
            
            % 고정될 node 좌표 정하기.
            if (ely-round((nely+2)/2))^2+(elx-round((nelx+2)/2))^2 <= 0.2*(rad)^2
                FixNode(end+1,:) = [ely,elx]; % node 좌표
                Fix(ely,elx) = 1; % binary 로 표시
                FixNumber(end+1) = (elx-1)*(nely+1)+ely; % node 넘버링
            end
            
            % 여기는 element 단위
            
            if ely ~= nely+1 && elx ~= nelx+1
                % Passive 영역 - 휠 바깥쪽 영역.
                if (ely-0.5-(nely)/2)^2+(elx-0.5-(nelx)/2)^2 > ((nely/2-0.5)-indent)^2 
                    Passive(ely,elx) = 1;
                else
                    Passive(ely,elx) = 0;
                end

                % Rim 부분은 항상 재료가 채워지도록 설정.
                %if ((ely-0.5-(nely)/2)^2+(elx-0.5-(nelx)/2)^2 <= (nely/2-0.5-indent)^2 && (ely-0.5-(nely)/2)^2+(elx-0.5-(nelx)/2)^2 > (floor(0.85*nely/2)-0.5-indent)^2)
                %    Active_Rim(ely,elx) = 1;
                %else
                %    Active_Rim(ely,elx) = 0;
                %end
            end
    end
end

piece=5
degree = 360/piece
rad = (degree*pi)/180
r = 64
% piece 부분 
for i =[1:64]
    for j = [64:128]
        if ((65-i)*(65-i)+(64-j)*(64-j)) < r^2
            if ((65-i)/(j-64)) < tan(rad)
                Passive(i, j) = 0;
            else
                Passive(i,j) = 1;
         end
        end
    end
end
for i = [64:128]
    for j = [1:128]
        Passive(i,j)= 1;
    end
end

for i = [1:64+1]
    for j =[1:64+1]
        Passive(i,j)=1;
    end
end
%% Load 
F = sparse(2*(nely+1)*(nelx+1),1); 
div_n = 90;
LoadSet = cell(div_n,1);
tempSet = cell(div_n,1);

for i = 1 : size(LoadBC,1)
  center = [round((nely+2)/2), round((nelx+2)/2)];
  tmp = LoadBC(i,:);
  n = (nely+1)*(tmp(2)-1)+tmp(1); 
  vec1(i,:) = [-(center(1)-tmp(1))/sqrt((center(1)-tmp(1))^2+(center(2)-tmp(2))^2) (center(2)-tmp(2))/sqrt((center(1)-tmp(1))^2+(center(2)-tmp(2))^2)];
  if tmp(1) < (nely+2)/2
      vec2(i,:) = null(vec1(i,:));
  else
      vec2(i,:) = -null(vec1(i,:));
  end


  for j = 1 : div_n
       if acos(dot(-vec1(1,:),-vec1(i,:)))*180/pi >= 180/div_n*(j-1) && acos(dot(-vec1(1,:),-vec1(i,:)))*180/pi < 180/div_n*j
           LoadSet{j}(end+1,:) = vec1(i,:);
           tempSet{j}(end+1,:) = vec2(i,:);
       elseif acos(dot(-vec1(1,:),-vec1(i,:)))*180/pi == 180/div_n*j
           LoadSet{j}(end+1,:) = vec1(i,:);
       end
  end

end


for i = 1 : size(LoadBC,1)
  tmp = LoadBC(i,:);
  n = (nely+1)*(tmp(2)-1)+tmp(1); 

  for j = 1 : div_n
       if acos(dot(-vec1(1,:),-vec1(i,:)))*180/pi >= 180/div_n*(j-1) && acos(dot(-vec1(1,:),-vec1(i,:)))*180/pi < 180/div_n*j
               F([2*n; 2*n-1],1) = 20*vec1(i,:)/size(LoadSet{j},1)'+w_load*20*vec2(i,:)/size(LoadSet{j},1)';
        elseif acos(dot(-vec1(1,:),-vec1(i,:)))*180/pi == 180/div_n*j
               F([2*n; 2*n-1],1) = 20*vec1(i,:)/size(LoadSet{j},1)'+w_load*20*vec2(i,:)/size(LoadSet{j},1)';
       end
  end
  
end


xTilde = x;
xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
loopbeta = 0;
cOld = inf;
c = 0;
%% START ITERATION
while loop < 100 && abs(cOld-c)>0.0001

  cOld = c;
  loopbeta = loopbeta + 1;
  loop = loop + 1;
%% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U = zeros(2*(nely+1)*(nelx+1),1);
  fixeddofs   = sortrows([FixNumber'*2;FixNumber'*2-1])';
  alldofs     = [1:2*(nely+1)*(nelx+1)];
  freedofs    = setdiff(alldofs,fixeddofs);
  U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
  U(fixeddofs,:)= 0;     
  
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  comp_hist(loop) = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)) + w*norm(Active-xPhys,1);
  for aa = 1 : size(xPhys,1)
      for bb = 1 : size(xPhys,2)
          if Active(aa,bb) > xPhys(aa,bb)
              dc_L1(aa,bb) = -1;
          else
              dc_L1(aa,bb) = 0;
          end          
      end
  end
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce + w*dc_L1;
  dv = ones(nely,nelx);
%% Compliance 에 대한 Sensitivity
  dx = beta*exp(-beta*xTilde)+exp(-beta);
  dc(:) = H*(dc(:).*dx(:)./Hs);
  dc = dc./max(max(abs(dc)));
  dv(:) = H*(dv(:).*dx(:)./Hs);
  dv = dv./max(max(abs(dv)));
%% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  l1 = 0; l2 = 1e+35; move = 0.15;
  while (l2-l1 > 1e-4)
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1.,min(x+move,x.*sqrt(abs(-dc./dv/lmid))))));
        
        %% 1차 density filter
        xTilde(:) = (H*xnew(:))./Hs;
        %% 2차 Threshold projection filter
        xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);  
        xPhys(find(Passive)) = 0;
        %xPhys(find(Active_Rim)) = 1;
          
        if sum(sum(xPhys)) - volfrac*nelx*nely > 0 % 
            l1 = lmid;
        else
            l2 = lmid;
        end      
      

  end  
  
  x_rs = reshape(xPhys,nelx*nely,1);
  Mnd2=sum(4*x_rs(:).*(1-x_rs(:)))/length(x_rs(:))*100; % Gray scale.
  change = max(max(abs(xnew-x)));
  x = xnew;
%  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
%       ' Vol.: ' sprintf('%6.3f',sum(sum(xPhys))/(nelx*nely)) ...
%        ' ch.: ' sprintf('%6.3f',change )])

% PLOT DENSITIES  
  if loop>=2

     if abs(Mnd2-Mnd2_old)<1 && beta<betaMax && abs(sum(sum(xPhys))/(nelx*nely)-volfrac) < 0.01
        beta=beta*1.5;
        loopbeta = 0;
     %   disp(['beta parameter is changed to ',num2str(beta)])    
   
    end
  end
  
  Mnd2_old=Mnd2;
  
 %% PRINT RESULTS
    
  colormap(gray); imagesc(-xPhys); axis equal; axis tight; axis off;pause(1e-6);
  if savefig
      filename = [pwd,filesep,sprintf('Iteration_%g',loop)]
      saveas(gcf,[filename,'.jpg']);
  end
   

end 


