function [S,D,xvalue]=GSE(infile,Lread,K,error)
if(nargin==3)
    err=0.005;
elseif(nargin==4)
    err=error;
else
    disp('Usage:GSE(infile,L,K,[error])');
    return;
end
close all
%figure;

M=load(infile);
x=M(:,1);
y=M(:,2);
%number of kmer 
s=sum(x.*y);

disp(['calculate the half peak width ...']);
ss=0;
for i=1:length(x);
    ss=ss+x(i)*y(i);
    if(ss/s>=err)
        L=i;
        break;
    end;
end
disp(['half peak width=' num2str(L)]);
disp('start peak detecting ...');
loc=[];
maxpeak=0;
maxloc=0;
for i=L+1:length(y)-L   
    ispeak=1;
    for j=i-L:i+L
        if(y(i)<y(j)||y(i)<mean(y)/5)
            ispeak=0;
            break;
        end
    end
    if(ispeak==1)
        loc=[loc i];
        if(y(i)>=maxpeak)
            maxpeak=y(i);
            maxloc=i;
        end
    end
end

minx=x(find(y(1:maxloc)==min(y(1:maxloc)),1));
xvalue=x(maxloc)+minx*minx/x(maxloc);
D=xvalue*Lread/(Lread-K+1);

S=s/(xvalue);

plot(x,y,'r');
hold on
line([xvalue xvalue],[0 1.3*max(y)]);

scatter(loc,y(loc),'r*')
axis([-1 28*x(maxloc) -1 12*y(maxloc)]);
title(infile);

% 
disp(['find ' num2str(length(loc)) ' peaks']);
disp(['the highest peak is (' num2str(x(maxloc)) ',' num2str(y(maxloc)) ')']);
disp(['estimate depth of original data is ' num2str(D)]);
disp(['estimate genome size of original data is ' num2str(S)]);
end
%    

