

for i=773:825
    format='[process check   (%d)   ]:\t-->\tok';
    format2='send to::%d.%d.%d.%d.%d';
    str=sprintf(format,i);
    str2=sprintf(format2,randi(500),randi(500),randi(500),randi(500),randi(500));
    disp(str2);
disp(str);
end
disp('connecting...');
disp('[host]--pralell:[tyama_lab:linux_os/tower 1 to 6]');
   disp('success!');