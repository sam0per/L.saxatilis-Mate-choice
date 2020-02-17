%This is the code for simulating the AS, SimMR, SimOR and SimM0 models in
%the article entitled "Assortative mating, sexual selection and their
%consequencs for the barrier to gene flow in Littorina", by S. Perini, M.
%Rafajlovic, A. M. Westram, K. Johannesson, and R. K. Butlin.
 

%Some model details are given next.
%In total, L loci are simulated.
%The divergently selected trait is determined by the first L/2 loci in
%females, whereas it is determined by the next L/2 loci in males.
%The phenotype is equal to the sum of allele-effect sizes at the
%corresponding L/2 loci in females/males.

%The initial condition is obtained from separate simulations with random
%mating, instead of assortative mating. The input files for the initial
%condition have the word "NULL" in their names (please see below). In
%this case, we had 200 such initial conditions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);
clear s

%%%%%%%%%%%%%%%%%
%Parameters.
%%%%%%%%%%%%%%%%%
%Set the ID of the simulation (ID is called partII here, and it is set to 1
%for simplicity, but it can take values from 1 to 200).

partII=1;

%Number of demes: M
M=400;
%Number of individuals in each deme
N=100;

%Mutatioon rate: mu (here, no mutations)
mu=0;
%Standard deviation for mutation-effect sizes: sigmaMuSel (here, set to 0,
%as there are no mutations)
sigmaMuSel=0;
%number of selected loci
L=80;
%Allow for a possibility that the phenotype of males is equal to othe sum
%of allele-effect sizes + some shift (equal to log(factor) here); but in
%this case, we set tthe shift to 0, by setting factor=1.
factor=1;

part=partII;
%Standard deviation for the environmental contribution to the phenotype: 
%sigmaE (but here, this is neglected, and sigmaE is set to 0).
sigmaE=0;
%Standard deviation for dispersal: sigmaD
sigmaD=1.5;
%Modelling Gaussian dispersal while accounting for the discreetness of the
%habitat:
x=-M:M;
pr=(erf((1-2*x)/(2*sqrt(2)*sigmaD))+erf((1+2*x)/(2*sqrt(2)*sigmaD)))/2;
pos0=kron((1:M)',ones(1,N/2));

%All loci are assumed to be unlinked.
r1=0.5;
r0=r1;
%neutr_time: set to 0 here, because selection acts from the start
neutr_time=0;
%Total running time with the assortment: sel_time
sel_time=10000;
%Total running time of the simulation: repeat
repeat=neutr_time+sel_time;

%Statistics are measured in intervals of "span" generations, here set to
%1000.
span=1000;

%Selection is divergent. Environmental transition occurs at deme denoted
%by: transition (here chosen deme in the middle of the habitat, i.e. deme
%M/2).
transition=M/2;


%Trait optima left (theta1) and right (theta2) of the environmental
%transition
theta1=2;
theta2=-theta1;
%Possible allele effect sizes at each locus are denoted by AlSize1 and
%AlSize2=-AlSize1. They are assumed to be the same for all loci, and they
%are chosen so that overshooting is possible.
%reached
AlSize1=2*theta1/(L);
AlSize2=2*theta2/(L);

%Total fitness disadvantage of an individual perfectly adapted to one
%habitat end when moved into the oother habitat end: sTot.
sTot=0.3;

%Translate sTot to the parameter sigma (inverse of the strengt of
%selection) in Fisher's geometric model 
sigma=sqrt(-2*theta1^2/log(1-sTot));





%Mating parameters: 
%Number of offspring per female: numMate
numMate=100;
count=0;
partI=part;
%Quantities of interest
clineAlF=zeros(repeat/span,M,L);%frequency of allele "AlSize1"
clineAlM=zeros(repeat/span,M,L);
clineAlAll=zeros(repeat/span,M,L);

LDF=zeros(repeat/span,M,L*(L-1)/2);
LDM=zeros(repeat/span,M,L*(L-1)/2);
LDALL=zeros(repeat/span,M,L*(L-1)/2);
LDFrel=zeros(repeat/span,M,L*(L-1)/2);
LDMrel=zeros(repeat/span,M,L*(L-1)/2);
LDALLrel=zeros(repeat/span,M,L*(L-1)/2);

PhenF=zeros(repeat/span,M,N/2);
PhenM=zeros(repeat/span,M,N/2);
pop1MALL=zeros(repeat/span,M,2*L,N/2);
pop1FALL=zeros(repeat/span,M,2*L,N/2);
AllMalesPhen=zeros(repeat/span,M);
matedMalesPhen=zeros(repeat/span,M);

%Starting condition: all homozygotes for AlSize1 in the first half of the
%habitat, and all homozygotes for AlSize2 in the second hald of the habitat

    pop1M=AlSize1*ones(M,2*L,N/2);
    pop1M(transition+1:M,:,:)=AlSize2;
    
    pop1F=AlSize1*ones(M,2*L,N/2);
    pop1F(transition+1:M,:,:)=AlSize2;
    Phen1F=reshape((sum(pop1F(:,1:L/2,:),2)+sum(pop1F(:,L+1:L+L/2,:),2)),M,[])+randn(M,N/2).*sigmaE;
    Phen1M=reshape((sum(pop1M(:,L/2+1:L,:),2)+sum(pop1M(:,L+L/2+1:2*L,:),2)),M,[])+log(factor)+randn(M,N/2).*sigmaE;
    
for I=1:repeat

    prog_num=zeros(M,N);
    
        sigmaS=sigma;
        sigmaMu=sigmaMuSel;
    
 
    pop_new1F=zeros(M,2*L,2*N);
    pop_new1M=zeros(M,2*L,2*N);
    
    Phen1F_new=zeros(M,2*N);
    Phen1M_new=zeros(M,2*N);

    num_popF=zeros(1,M);
    num_popM=zeros(1,M);


    %migration of males
    a=rand(M,N/2);
    [~,b]=histc(a,cumsum([0 pr]));
    pos1=pos0+b-M-1;
    %now assure that those that would overshoot the habitat size actually
    %go to patch 1 or M; depending on where the overshooting occurs
    a1=find(pos1<1);
    pos1(a1)=1;
    a1=find(pos1>M);
    pos1(a1)=M;
    
   
    
   
    
    for j=1:M
            a1=find(pos1==j);
            num_popM(1,j)=length(a1);
            k=ceil(a1/M);
            i=a1-floor(a1/M)*M;
            as=find(i==0);
            i(as)=M;
            for i1=1:num_popM(1,j)
                pop_new1M(j,:,i1)=reshape(pop1M(i(i1),:,k(i1)),[],1);
                Phen1M_new(j,i1)=Phen1M(i(i1),k(i1));

            end
    end
   
    
      %migration of females
    a=rand(M,N/2);
    [~,b]=histc(a,cumsum([0 pr]));
    pos1=pos0+b-M-1;
    %now assure that those that would overshoot the habitat size actually
    %go to patch 1 or M; depending on where the overshooting occurs
    a1=find(pos1<1);
    pos1(a1)=1;
    a1=find(pos1>M);
    pos1(a1)=M;
   
    
    for j=1:M
            a1=find(pos1==j);
            num_popF(1,j)=length(a1);
            k=ceil(a1/M);
            i=a1-floor(a1/M)*M;
            as=find(i==0);
            i(as)=M;
            for i1=1:num_popF(1,j)
                pop_new1F(j,:,i1)=reshape(pop1F(i(i1),:,k(i1)),[],1);
                Phen1F_new(j,i1)=Phen1F(i(i1),k(i1));
            end
    end
    
    
    
 
   

    %now ask who mates with whom
    %allow for finite number of mating attempts
    
    AllMalesPhen1=[];
    for i=1:M
        %all mating pairs possible
        a=zeros(1,num_popF(1,i)*num_popM(1,i));
        b=Phen1F_new(i,1:num_popF(1,i));
        b=reshape(b,1,[])';
        b=kron(b,ones(1,num_popM(1,i)));
        c=Phen1M_new(i,1:num_popM(1,i));
        c=kron(c,ones(num_popF(1,i),1));
        ratio=b-c;
pI=ones(length(ratio(:,1)),length(ratio(1,:)));        
%now compute the total propensity of each female to mate by summing
        %pI over all available males in the habitat
       pI=pI';
        a=sum(pI);
        a=kron(a,ones(length(pI(:,1)),1));
        pI=pI./a;
        pI(end,:)=2;
boxesMate=cumsum([zeros(1,length(pI(1,:)));pI]);
        a1=(1:num_popF(1,i))*10;
        a1=kron(a1,ones(length(boxesMate(:,1)),1));
        boxesMate=boxesMate+a1;%a1=reshape(a1,1,[]);
        boxesMate=reshape(boxesMate,1,[]);
                %now mating
        a=rand(1,num_popF(1,i)*numMate);%this is to assign probabilities for each mating
        
        a2=(1:num_popF(1,i))*10;
        a2=kron(a2,ones(numMate,1));
        a=a+reshape(a2,1,[]);
        %now ask where each male falls in
        [~,b]=histc(a,boxesMate);        b1=(0:(num_popF(1,i)-1))*(num_popM(1,i)+1);
        b1=kron(b1,ones(numMate,1));
        b1=reshape(b1,1,[]);
        b=b-b1;
        %this is now the list of males
        males=reshape(pop_new1M(i,:,1:num_popM(1,i)),2*L,[]);
        
        malesMate=males(1:L,b);
        matedMalesPhen1(1,i)=mean(Phen1M_new(i,b),2);
        AllMalesPhen1(1,i)=mean(Phen1M_new(i,1:num_popM(1,i)),2);

        malesMate1=males(L+1:2*L,b);
        %choose alleles that they contribute with
        a=rand(L,length(b));
        cc=[0 0.5 1.1];
        [~,cc1]=histc(a,cc);
        cc1=find(cc1==2);
        malesMate(cc1)=malesMate1(cc1);
        
        %now for females
        a=kron((1:num_popF(1,i)),ones(numMate,1));
        a=reshape(a,1,[]);
        females=reshape(pop_new1F(i,:,1:num_popF(1,i)),2*L,[]);
        femalesMate=females(1:L,a);
        femalesMate1=females(L+1:2*L,a);
        a=rand(L,length(a));
        cc=[0 0.5 1.1];
        [~,cc1]=histc(a,cc);
        cc1=find(cc1==2);
        femalesMate(cc1)=femalesMate1(cc1);
        prog1=[femalesMate;malesMate];

        %now apply soft viability selection to keep the number of offspring
        %the same everywhere
        %here you can put in mutations if necessary
        
        
        %compute fitnesses
 %but since selection acts on phenotypes, that are now determined
        %differently in males and females, must already here assign sex!
        %take a random permutation of offspring, and take the first half to
        %be females; the second half are all males
        a=randperm(length(prog1(1,:)));
        prog1=prog1(:,a);

        prog1F=prog1(:,1:floor(length(prog1(1,:))/2));
        prog1M=prog1(:,length(prog1F(1,:))+1:end);
        phFactorF=sum(prog1F(1:L/2,:),1)+sum(prog1F(L+1:L+L/2,:),1);%these are females
        phFactorM=(sum(prog1M(L/2+1:L,:),1)+sum(prog1M(L+L/2+1:2*L,:),1))+log(factor);%these are males
        
        
%start with females
        
        
        a=phFactorF;
        a=a+randn(length(a(:,1)),length(a(1,:))).*sigmaE;
        phen=a;
        if i<=transition
            fit1=exp(-(a-theta1).^2/(2*sigmaS^2));
        else
            fit1=exp(-(a-theta2).^2/(2*sigmaS^2));
        end
        
        a1=cumsum([0 fit1/sum(fit1)]);
        a1(end)=2*a1(end);
        %now choose N/2 females only
        a=rand(1,N/2);
        [~,b]=histc(a,a1);
        c1=prog1F(:,unique(b));
        phen1=phen(:,unique(b));
        phen2=phen;
        c2=prog1;
        c2(:,unique(b))=[];
        phen(:,unique(b))=[];
        fit1(unique(b))=[];
        if length(unique(b))~=N/2
            a=rand(1,N/2-length(unique(b)));
            a1=cumsum([0 fit1/sum(fit1)]);
            a1(end)=2*a1(end);
            [~,b]=histc(a,a1);
            c1=[c1 c2(:,b)];
            phen1=[phen1 phen(:,b)];
        end
        Phen1F(i,:)=phen1;
        pop1F(i,:,:)=c1;
        
        %now repeat for males
        a=phFactorM;
        a=a+randn(length(a(:,1)),length(a(1,:))).*sigmaE;
        phen=a;
        if i<=transition
            fit1=exp(-(a-theta1).^2/(2*sigmaS^2));
        else
            fit1=exp(-(a-theta2).^2/(2*sigmaS^2));
        end
        
        a1=cumsum([0 fit1/sum(fit1)]);
        a1(end)=2*a1(end);
        %now choose N/2 males only
        a=rand(1,N/2);
        [~,b]=histc(a,a1);
        c1=prog1M(:,unique(b));
        phen1=phen(:,unique(b));
        phen2=phen;
        c2=prog1;
        c2(:,unique(b))=[];
        phen(:,unique(b))=[];
        fit1(unique(b))=[];
        if length(unique(b))~=N/2
            a=rand(1,N/2-length(unique(b)));
            a1=cumsum([0 fit1/sum(fit1)]);
            a1(end)=2*a1(end);
            [~,b]=histc(a,a1);
            c1=[c1 c2(:,b)];
            phen1=[phen1 phen(:,b)];
        end
        Phen1M(i,:)=phen1;
        pop1M(i,:,:)=c1;
        
        
        
    end
    %nowdone with all parts in one lifecycle.
    %time to sample populations.
 
    
if mod(I,span)==0

count=count+1;

PhenF(count,:,:)=Phen1F; 
pop1FALL(count,:,:,:)=pop1F;
pop1MALL(count,:,:,:)=pop1M;


PhenM(count,:,:)=Phen1M;          
a=zeros(M,L,N);
a(:,:,1:N/2)=pop1F(:,1:L,:);
a(:,:,N/2+1:N)=pop1F(:,L+1:2*L,:);
a1=a;
b=find(a(:,:,:)<0);
a(b)=0;
b=(a(:,:,:)-floor(a)==0);
b=find(b==0);
a(b)=1;
a=sum(a,3)/length(a(1,1,:));
a=reshape(a,M,[]);
clineAlF(count,:,:)=a;



a(:,:,1:N/2)=pop1M(:,1:L,:);
a(:,:,N/2+1:N)=pop1M(:,L+1:2*L,:);
a2=a;
b=find(a(:,:,:)<0);
a(b)=0;
b=(a(:,:,:)-floor(a)==0);
b=find(b==0);
a(b)=1;
a=sum(a,3)/length(a(1,1,:));
a=reshape(a,M,[]);
clineAlM(count,:,:)=a;  

b1=zeros(M,L,2*N);
b1(:,:,1:N)=a1;
b1(:,:,N+1:2*N)=a2;
a=b1;
b=find(a(:,:,:)<0);
a(b)=0;
b=(a(:,:,:)-floor(a)==0);
b=find(b==0);
a(b)=1;
a=sum(a,3)/length(a(1,1,:));
a=reshape(a,M,[]);
clineAlAll(count,:,:)=a;
for i=1:M
    cc1=0;

    for j=1:L-1
        for k=(j+1):L
            cc1=cc1+1;
            a=[[reshape(pop1F(i,j,:),1,[]) reshape(pop1F(i,L+j,:),1,[])];[reshape(pop1F(i,k,:),1,[]) reshape(pop1F(i,k+L,:),1,[])]];
            b=[[reshape(pop1M(i,j,:),1,[]) reshape(pop1M(i,L+j,:),1,[])];[reshape(pop1M(i,k,:),1,[]) reshape(pop1M(i,k+L,:),1,[])]];
            a=ceil(a);
            b=ceil(b);
            c=[a b];
            
            a1=sum(a);
            a11=length(find(a1==2));
            a22=length(find(a1==0));
            if length(find(a1==1))>=1
                a12=length(find(a1==1 & a(1,:)==1))/length(a(1,:));
                a21=length(find(a1==1 & a(1,:)==0))/length(a(1,:));
            else
                a12=0;
                a21=0;
            end
            p1=length(find(a(1,:)==1))/length(a(1,:));
            p2=length(find(a(2,:)==1))/length(a(1,:));
            LDF(count,i,cc1)=a11*a22-a12*a21;
            if p1*(1-p1)*p2*(1-p2)~=0
                LDFrel(count,i,cc1)=LDF(count,i,cc1)/sqrt(p1*(1-p1)*p2*(1-p2));
            else
                LDFrel(count,i,cc1)=NaN;
            end
            
            a=b;
            a1=sum(a);
            a11=length(find(a1==2))/length(a(1,:));
            a22=length(find(a1==0))/length(a(1,:));
            if length(find(a1==1))>=1
                a12=length(find(a1==1 & a(1,:)==1))/length(a(1,:));
                a21=length(find(a1==1 & a(1,:)==0))/length(a(1,:));
            else
                a12=0;
                a21=0;
            end
            p1=length(find(a(1,:)==1))/length(a(1,:));
            p2=length(find(a(2,:)==1))/length(a(1,:));
            LDM(count,i,cc1)=a11*a22-a12*a21;
            if p1*(1-p1)*p2*(1-p2)~=0
                LDMrel(count,i,cc1)=LDM(count,i,cc1)/sqrt(p1*(1-p1)*p2*(1-p2));
            else
                LDMrel(count,i,cc1)=NaN;
            end
            
            a=c;
            a1=sum(a);
            a11=length(find(a1==2))/length(a(1,:));
            a22=length(find(a1==0))/length(a(1,:));
            if length(find(a1==1))>=1
                a12=length(find(a1==1 & a(1,:)==1))/length(a(1,:));
                a21=length(find(a1==1 & a(1,:)==0))/length(a(1,:));
            else
                a12=0;
                a21=0;
            end
            p1=length(find(a(1,:)==1))/length(a(1,:));
            p2=length(find(a(2,:)==1))/length(a(1,:));
            LDALL(count,i,cc1)=a11*a22-a12*a21;
            if p1*(1-p1)*p2*(1-p2)~=0
                LDALLrel(count,i,cc1)=LDALL(count,i,cc1)/sqrt(p1*(1-p1)*p2*(1-p2));
            else
                LDALLrel(count,i,cc1)=NaN;
            end
            
            
        end
    end
end
matedMalesPhen(count,:)=matedMalesPhen1;
AllMalesPhen(count,:)=AllMalesPhen1;

end
end


filename=sprintf('AssortMatingNULL_sTot%g_N%g_M%g_L%g_r%g_sigmaD%g_theta%g_numMate%g_sigmaE%g_factor%g_part%g', sTot, N, M, L, r1, sigmaD, theta1, numMate, sigmaE, factor, partI);

filename=strrep(filename,'.','p');
save(filename,'LDF','LDM','LDALL','LDFrel','LDMrel','LDALLrel','matedMalesPhen','AllMalesPhen','factor','pop1FALL','pop1MALL','transition','AlSize1','AlSize2','sTot','sigma','sel_time','sigmaMuSel','r1','r0','L','repeat','sigmaD','N','sigmaMu','mu','span','M','theta1','theta2', 'numMate','clineAlF','clineAlAll', 'clineAlM','PhenF','PhenM', 'sigmaE');            

%repeat this simulation (with different partII) to obtain an ensemble of
%independent initial conditions (that will be used as an input for the
%models with assortment).
  


