function [SMF,UR,UL,SMFsym,URsym,ULsym,del,Ps,D,eps] = smform(G)
%SMFORM  Convert transfer function matrix to Smith-McMillan Form
%
%   SMF = SMFORM(G) returns the equivalent Smith-McMillan form
%   matrix SMF of transfer function matrix G as a transfer functino object.
%   G is intended to be a transfer function object, but should also work as
%   a state-space or ZPK object.
% 
%   [SMF,UR,UL] = SMFORM(G) returns the SMF as well as left and right
%   unimodular transfer function matrics UL and UR, but only if G is a
%   square matrix. If G is not square then UR and UL are empty matrices.
%
%   [SMF,UR,UL,SMFsym,URsym,ULsym,del,Ps,D,eps] = SMFORM(G) returns the SM
%   form and unimodular matrices as both transfer function objects and as
%   symbolic objects. Additionally, computational steps for the SM form are
%   provided as symbolic objects:
%       - greatest common denominator 'del'
%       - polynomial matrix 'Ps'
%       - determinantal divisors 'D'
%       - epsilon terms 'eps'

%% documentation
%   author: Kevin R. Mallon
%   email: krmallon@ucdavis.edu
%   created: April 2015;
%
%   Recurring commands in code:
% 
%       var=prod(factor(var)) -- Factors the symbolic expression 'var' in
%           order to ensure cancelling of terms. Equivalent to minreal, but
%           for symbolic expressions.

%% begin
s=tf('s');
syms x
G=minreal(G);
[r,m]=size(G);

%% convert from transfer function to symbolic
[TFMnum,TFMden]=tfdata(G);
Gn=x*zeros(r,m);
Gd=x*zeros(r,m);
for i=1:r
    for j=1:m
        num=cell2mat(TFMnum(i,j));
        den=cell2mat(TFMden(i,j));
        Gn(i,j)=poly2sym(num,x);
        Gd(i,j)=poly2sym(den,x);
    end
end

%% determine denominator (del)
del=lcm(reshape(Gd,[1,numel(Gd)]));
del=prod(factor(del));
g_del=sym2poly(del);
del=del/g_del(1);

%% form polynomial matrix
Ps=del*(Gn./Gd);
for i=1:r
    for j=1:m
        Ps(i,j)=prod(factor(Ps(i,j)));
    end
end

%% determine D_k
D0=1;D=x*zeros(1,min(r,m));
for k=1:min(r,m)
    clear M
    % minor of order k
        rowstoclear=r-k;
        columnstoclear=m-k;

    % establish combinations rows/columns to cancel 
        allrows=1:r;
        allcolumns=1:m;
        rowcombs=allrows;
        colcombs=allcolumns;
        for i=1:rowstoclear-1
            rowcombs=combvec(rowcombs,allrows(i+1:end));
        end
        for i=1:columnstoclear-1
            colcombs=combvec(colcombs,allcolumns(i:end));
        end
        [~,numrcom]=size(rowcombs);
        [~,numccom]=size(colcombs);

    % remove invalid row combinations & sort
        i=1;
        while i<=numrcom
            W=rowcombs(:,i);
            if issorted(W)==0
                rowcombs(:,i)=[];
            elseif numel(unique(W))<numel(W)
                rowcombs(:,i)=[];
            else
                i=i+1;
            end
            [~,numrcom]=size(rowcombs);
        end
        rowcombs=sort(rowcombs,'descend');

    % remove invalid column combinations & sort
        i=1;
        while i<=numccom
            W=colcombs(:,i);
            if issorted(W)==0
                colcombs(:,i)=[];
            elseif numel(unique(W))<numel(W)
                colcombs(:,i)=[];
            else
                i=i+1;
            end
            [~,numccom]=size(colcombs);
        end
        colcombs=sort(colcombs,'descend');

    % evaluate minor determinants
        for i=1:numrcom
            for j=1:numccom
                Mw=Ps;
                RC=rowcombs(:,i);
                CC=colcombs(:,j);
                for p=1:rowstoclear
                    Mw(RC(p),:)=[];
                end
                for p=1:columnstoclear
                    Mw(:,CC(p))=[];
                end
                M(i,j)=det(Mw);
            end
        end
    
    % find greatest common denominator
        GCD_M=gcd(reshape(M,[1,numel(M)]));
        g_gcdm=sym2poly(GCD_M);
        GCD_M=GCD_M/g_gcdm(1);
        
    % record value for D_k    
        D(k)=prod(factor(GCD_M));
end
D=[D0,D];

%% determine epsilons
eps=x*zeros(1,length(D)-1);
for i=1:length(D)-1
    eps(i)=prod(factor(D(i+1)/D(i)));
end

%% determine symbolic Smith-McMillan Form matrix
SMFsym=x*zeros(r,m);
for i=1:length(eps);
    SMFsym(i,i)=prod(factor(eps(i)/del));
end

%% determine unimodular matrices
% matrix fraction decomposition
% G = inv(UL)*SMF*inv(UR)
if r==m
    Nr=Ps;
    Dr=del*eye(m);

    % isolate numerators, denominators
    for i=1:length(eps);
        SMFii=prod(factor(SMFsym(i,i)));
        [num,den]=numden(SMFii);
        smf_num(i)=prod(factor(num));
        smf_den(i)=prod(factor(den));
    end
    Nd=diag(smf_num);
    Dd=diag(smf_den);

    isunimodular=0;
    while isunimodular==0
        URsym=Dr*Dd^-1;
        if length(sym2poly(det(URsym))) == 1
            isunimodular=1;
        else
            Nr=Nr*URsym^-1;
            Dr=Dr*URsym^-1;
        end
    end
    % URsym=Dr*Dd^-1;
    ULsym=(Nr*Nd^-1)^-1;

    % isunimodular=0;
    % while isunimodular==0
    %     ULsym=(Nr*Nd^-1)^-1
    %     if length(sym2poly(det(ULsym))) == 1
    %         isunimodular=1;
    %     else
    %         Nr=ULsym^-1*Nr;
    %         Dr=ULsym^-1*Dr;
    %     end
    % end
    % URsym=Dr*Dd^-1;
else
    Nr=Ps;
    Dr=del*eye(m);

    % isolate numerators, denominators
    for i=1:length(eps);
        SMFii=prod(factor(SMFsym(i,i)));
        [num,den]=numden(SMFii);
        smf_num(i)=prod(factor(num));
        smf_den(i)=prod(factor(den));
    end
    Nd=diag(smf_num);
    Dd=diag(smf_den);

    isunimodular=0;
    while isunimodular==0
        URsym=Dr*Dd^-1;
        if length(sym2poly(det(URsym))) == 1
            isunimodular=1;
        else
            Nr=Nr/URsym;
            Dr=Dr/URsym;
        end
    end
    %URsym=Dr*Dd^-1;
    %ULsym=(Nr*Nd^-1)^-1;
end
%% convert from symbolic to transfer function
SMF=0.*tf(G);
for i=1:r
    for j=1:m
        [num,den]=numden(SMFsym(i,j));
        num=sym2poly(num);
        den=sym2poly(den);
        SMF(i,j)=tf(num,den);
    end
end

if r == m
    UR=minreal(zeros(size(URsym))*s/s);
    for i=1:length(UR)
        for j=1:length(UR)
            [num,den]=numden(URsym(i,j));
            num=sym2poly(num);
            den=sym2poly(den);
            UR(i,j)=tf(num,den);
        end
    end

    UL=minreal(zeros(size(ULsym))*s/s);
    for i=1:length(UL)
        for j=1:length(UL)
            [num,den]=numden(ULsym(i,j));
            num=sym2poly(num);
            den=sym2poly(den);
            UL(i,j)=tf(num,den);
        end
    end
else
    UR = [];
    UL = [];
end
%% determine poles and zeros
% [zz,pp,~] = zpkdata(SMF);
% Z=[];P=[];
% for i=1:length(SMF)
%     Z=[Z;cell2mat(zz(i,i))];
%     P=[P;cell2mat(pp(i,i))];
% end
% Z=sort(Z);
% P=sort(P);

%% reformat symbolic expressions
syms s
SMFsym=subs(SMFsym,x,s);
del=subs(del,x,s);
Ps=subs(Ps,x,s);
D=subs(D,x,s);
eps=subs(eps,x,s);

if r==m
    URsym=subs(URsym,x,s);
    ULsym=subs(ULsym,x,s);
end

end


