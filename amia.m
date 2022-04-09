
clear
supsize=10;
memsize=5;
size=500;
supcnt=1;
gene_size=100;
stops=16;
p_elt=0.05;
x=1;
y=2;
seed=0;
template=zeros(1,stops);
simsup=zeros(1,stops);
sim=zeros(size,stops);
city=zeros(stops,2);
fitness=zeros(1,gene_size);
plt_fit=zeros(1,gene_size);
elite=zeros(p_elt*size,stops);%各世代のエリート
c=zeros(size,1);%解の濃度
ay=zeros(size,size);%解の共通度合い
ac=zeros(size,size);%TAC1を越えた個体
as=zeros(supsize,1);
asy=zeros(supsize,1);
normalize_e=zeros(size,1);%正規化した残りやすさ
TC=zeros(size,1);
TAC1=zeros(size,1);
TAC2=zeros(size,1);
TC_power=1.5;
e=ones(size,1);%解の残りやすさ
supress=zeros(supsize,stops);
memory=zeros(memsize,stops);
MemoryT=0.8;
lc_flg=false;
for i=1:size
    TC(i)=0.5;
    TAC1(i)=0.9;
end


%抗体の実体化
rng(seed,'twister');
for i=1:size
antb(i)=antibody(randperm(stops),'normal');
end


%都市集合の生成
for cnt = 1:stops
    theta=cnt*(2*pi)/stops;
    city(cnt,:)=[cos(theta),sin(theta)];
end

for g=1:gene_size
           template=zeros(1,stops);
           sim=zeros(size,stops);
   
    seed=seed+1;
rng(seed,'twister');
    normalize_e=e/sum(e);
 %step2:交換交叉
 flag=true;
for i=1:size*0.1
while flag==true
    rs=find(normalize_e<rand);
    if numel(rs)~=0  
        slc1=rs(randi(numel(rs)));
    else
        slc1=randi(size);
    end   
    rs=find(normalize_e<rand);
    if numel(rs)~=0
        slc2=rs(randi(numel(rs)));
    else
         slc2=randi(size);
    end
    
   pos1=randi(stops);
   pos2=randi([pos1 stops]);
   tmp=antb(slc1).data;
   
   set1=antb(slc1).data(pos1:pos2);
   set2=antb(slc2).data(pos1:pos2);

   if(isempty(setdiff(set1,set2))==true)%差集合が0であれば入れ替える
        for i=pos1:pos2
            antb(slc1).data(i)=antb(slc2).data(i);
            antb(slc2).data(i)=tmp(i);
        end
        flag=false;
   end
end

%step3:突然変異
   seed=seed+1;
   rng(seed,'twister');
if randi(10000)<=100
   pos1=randi(stops);
   pos2=randi(stops);
           seed=seed+1;
        rng(seed,'twister');
    rs=find(normalize_e<=rand);
    if numel(rs)~=0  
        slc=rs(randi(numel(rs)));
    else
        slc=randi(size);
    end   
   tmp=antb(slc).data(pos1);
   antb(slc).data(pos1)=antb(slc).data(pos2);
   antb(slc).data(pos2)=tmp;
end
end

for i=1:size
dist=0;
for root=2:stops
    dist=dist+norm(city(antb(i).data(root),:)-city(antb(i).data(root-1),:));
end
dist=dist+norm( city(antb(i).data(stops),:)-city(antb(i).data(1),:));
fitness(i)=1./dist;
end
[fitmx,mx_i]=maxk(fitness(:),p_elt*size);
plt_fit(g)=fitmx(1);


for i=1:p_elt*size
   elite(i,:)=antb(mx_i(i)).data;
end
%disp('[エリート保存完了]');
n=mx_i(1);


%濃度計算
for v=1:size
    for w=1:size
        cmn_path=0;
       if v~=w
        cmn_path=stops-nnz(antb(v).data-antb(w).data);
        ay(v,w)=cmn_path./stops;
        if(ay(v,w)>TAC1)
            ac(v,w)=1;
        else
            ac(v,w)=0;
        end
       end
    end
end

for v=1:size
 c(v)=sum(ac(v,:))/size;
end

[hoge,I]=maxk(c,1);
src=1;
supcnt=1;
for v=1:size
    if c(v)>TC(v)
        if TC(v)>0.5*TC_power
          disp("sup");
           TC(v)=0.5;
           TAC1(v)=0.9;
           antb(v).type="sup";
           %抑制細胞群を作成
           if supcnt<=supsize
            supress(supcnt,:)=antb(I).data;%最も濃度の高い個体をサプレッサー細胞に保存
            supcnt=supcnt+1;
           else
           for k=1:supsize  
                cmn_path=0;
                for s=1:stops
                    if supress(k,s)==antb(I).data(s)
                        cmn_path=cmn_path+1;
                    end
                end
                asy(k)=cmn_path/stops;
           end
         [hoge,fuga]=maxk(as,1);
          supress(fuga,:)=antb(I).data;
           end 
         %抑制細胞作成完了 
        else
          disp('mem');
           %記憶細胞の作成
         antb(v).type="mem";
         TC(v)=TC(v)+normalize_e(v);
         TAC1(v)=TAC1(v)+normalize_e(v);
        end
    end
end


for v=1:supsize
    if asy(v)>TAC2
        as(v)=1;
    else
        as(v)=0;
    end
end

for v=1:size
    if ay(I,v)>TAC1(v)
            sim(src,:)=antb(I).data;
            src=src+1;
    end   
end


for s=1:stops
    if sum(sim([1:src-1],s))/(src-1)==sim(1,s)
        template(s)=sim(1,s);
    else
        template(s)=0;
    end
end

if sum(template)~=0
for mt=1:stops-fix(stops*MemoryT)
    seed=seed+1;
    rng(seed,'twister');
    template(randi(stops))=0;
end
end

%template
%{
for s=1:stops
    if sum(supress(:,s))/supsize==supress(1,s)
        simsup(s)=supress(1,s);
    else
        simsup(s)=0;
    end
end
%}

%次世代の作成
for v=1:size
    if strcmp(antb(v).type,"mem")==true
        %disp('mem');
        %記憶細胞によってつくられる抗体
         antb(v).data=template;
        for s=1:stops
            if antb(v).data(s)==0 
                r=setdiff([1:stops],antb(v).data);
                seed=seed+1;
                rng(seed,'twister');
                antb(v).data(s)=r(randi(numel(r)));
            end
        end
    elseif strcmp(antb(v).type,"sup")==true%抑制
        %disp('sup');       
        seed=seed+1;
        rng(seed,'twister');
        for ss=1:size
              antb(ss).data=randperm(stops);
              TAC2(ss)=TAC2(ss)-normalize_e(ss);
        end
    end
end


for v=1:size
    e(v)=fitness(v)*(prod(1-as));
end


for i=1:p_elt*size 
    antb(randi(size)).data=elite(i,:);
end


%経路の表示---------------------------
clf reset

cnt=1:stops;
 plot(city(cnt,x),city(cnt,y),'*','MarkerSize',5,'MarkerEdgeColor','black', 'MarkerFaceColor','black')
 axis equal
hold on
for cnt=1:stops-1
   plot([city(antb(n).data(cnt),x),city(antb(n).data(cnt+1),x)],[city(antb(n).data(cnt),y),city(antb(n).data(cnt+1),y)],'r','color','black');
   plot([city(antb(n).data(stops),x),city(antb(n).data(1),x)],[city(antb(n).data(stops),y),city(antb(n).data(1),y)],'r','color','black');
   drawnow
end
%------------------------------------
text=sprintf('[generation]::%d',g);
disp(text);
end%世代の終了

antb(n).data
%都市の表示
figure(1); 
cnt=1:stops;
    plot(city(cnt,x),city(cnt,y),'*','MarkerSize',5,'MarkerEdgeColor','black', 'MarkerFaceColor','black')
    axis equal
    hold on;
%経路の表示
cnt=1:stops-1;
   plot([city(antb(n).data(cnt),x),city(antb(n).data(cnt+1),x)],[city(antb(n).data(cnt),y),city(antb(n).data(cnt+1),y)],'r','color','black');
   plot([city(antb(n).data(stops),x),city(antb(n).data(1),x)],[city(antb(n).data(stops),y),city(antb(n).data(1),y)],'r','color','black');
   hold off
%{
figure(2);
i=1:gene_size;
    plot(i,plt_fit(i));
%}
xlabel('generation');
ylabel('fitness');

xp=1:gene_size;
ansplt=[xp; plt_fit(xp)];
fileID = fopen('amia.txt','w');
fprintf(fileID,'%f\t%f\n',ansplt);
fclose(fileID);
