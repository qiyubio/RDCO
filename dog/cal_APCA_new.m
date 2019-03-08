  
% Author: Qi Yu
% Contact: dryuqi@gmail.com
function cal_APCA_new(filename)
    
    M1n = dlmread(filename);

    %logm=log(M1n);
    %logm=log10(M1n - min(M1n));
    %logm=zscore(M1n);
    %logm=normc(M1n);
    
    logm=log2((M1n - min(M1n)) / ( max(M1n) - min(M1n)));
    %logm=log((M1n - min(M1n))/ ( max(M1n) - min(M1n)));
    %[logm,f]=log2(M1n  - min(M1n));
    
    %figure8=figure;
    %hold all;
    %plot(M1n);
    
    %figure9=figure;
    %plot(logm);
    %legend([{'Original data'} {'log2 transform'}]);
    %hold off;
    
    M1n = logm;
    fname=sprintf('%s_apca',filename);
    [retMatrix1,retMatrix2,retMatrix3,retMatrix4] = WithAPCA(M1n,30,20,fname);
    m=mean(logm(~isinf(logm)));
    
    c = zeros(size(retMatrix1.apca));
    [a,b]=size(retMatrix1.apca);
    
    
    c2 = zeros(size(retMatrix2.apca));
    [a2,b2]=size(retMatrix2.apca);
    
    
    c3 = zeros(size(retMatrix3.apca));
    [a3,b3]=size(retMatrix3.apca);
    
    
    c4 = zeros(size(retMatrix4.apca));
    [a4,b4]=size(retMatrix4.apca);
    
    
    for i = 1:b
        new=retMatrix1.apca(1,i).y-m;
        
        c(1,i)=new;
        
    end
    new=[retMatrix1.apca.lx;retMatrix1.apca.rx;retMatrix1.apca.y;c];
    
    fname=sprintf('%s_apca_best1_%d.file',filename,b);
    fid = fopen(fname, 'w');
    fprintf(fid,'%d\t%d\t%d\t%d\n',new);
    fclose(fid);
    
    for i = 1:b2
        new=retMatrix2.apca(1,i).y-m;
        
        c2(1,i)=new;
        
    end
    new=[retMatrix2.apca.lx;retMatrix2.apca.rx;retMatrix2.apca.y;c2];
    
    
    fname=sprintf('%s_apca_best2_%d.file',filename,b2);
    fid = fopen(fname, 'w');
    fprintf(fid,'%d\t%d\t%d\t%d\n',new);
    fclose(fid);
   
    for i = 1:b3
        new=retMatrix3.apca(1,i).y-m;
        
        c3(1,i)=new;
        
    end
    new=[retMatrix3.apca.lx;retMatrix3.apca.rx;retMatrix3.apca.y;c3];
    
    
    fname=sprintf('%s_apca_best3_%d.file',filename,b3);
    fid = fopen(fname, 'w');
    fprintf(fid,'%d\t%d\t%d\t%d\n',new);
    fclose(fid);
    
     for i = 1:b4
        new=retMatrix4.apca(1,i).y-m;
        
        c4(1,i)=new;
        
    end
    new=[retMatrix4.apca.lx;retMatrix4.apca.rx;retMatrix4.apca.y;c4];
    
    
    fname=sprintf('%s_apca_best4_%d.file',filename,b4);
    fid = fopen(fname, 'w');
    fprintf(fid,'%d\t%d\t%d\t%d\n',new);
    fclose(fid);
end

