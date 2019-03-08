function segment = apca(data,num_segments,dummy)

left_x = [1 : 2 : size(data,1)-1];
right_x      = left_x + 1;
right_x(end) = size(data,1);
number_of_segments = length(left_x );

for i = 1 : number_of_segments 
   
   segment(i).lx = left_x(i);
   segment(i).rx = right_x(i);
   segment(i).mc = inf;
   segment(i).y = inf;
end;

for i = 1 : number_of_segments - 1
   coef = polyfit([segment(i).lx :segment(i+1).rx]',data(segment(i).lx :segment(i+1).rx),0);
   best = (0*( [segment(i).lx :segment(i+1).rx]' ))+coef(1);
   segment(i).mc = sum((data([segment(i).lx :segment(i+1).rx]')-best).^2);
end;

while  length(segment) > num_segments  
      
   [value, i ] = min([segment(:).mc]);
      
     if i > 1 && i < length(segment) -1								
      
    	coef = polyfit([segment(i).lx :segment(i+2).rx]',data(segment(i).lx :segment(i+2).rx),0);        
    	best = (0 *( [segment(i).lx :segment(i+2).rx]' ))+coef(1);       
        segment(i).mc = sum((data([segment(i).lx :segment(i+2).rx]')-best).^2);
	 	segment(i).rx = segment(i+1).rx;
    	segment(i+1) = [];
       	i = i - 1;    
        coef = polyfit([segment(i).lx :segment(i+1).rx]',data(segment(i).lx :segment(i+1).rx),0);
      	best = (0*( [segment(i).lx :segment(i+1).rx]' ))+coef(1);    
        segment(i).mc = sum((data([segment(i).lx :segment(i+1).rx]')-best).^2);
       
   elseif i == 1
       
	   coef = polyfit([segment(i).lx :segment(i+2).rx]',data(segment(i).lx :segment(i+2).rx),0);
       best = (0*( [segment(i).lx :segment(i+2).rx]' ))+coef(1);
       segment(i).mc = sum((data([segment(i).lx :segment(i+2).rx]')-best).^2);
	   segment(i).rx = segment(i+1).rx;
      segment(i+1) = [];
              
   else
     
      segment(i).rx = segment(i+1).rx;
      segment(i).mc = inf;
      segment(i+1) = [];    
      i = i - 1;       
      coef = polyfit([segment(i).lx :segment(i+1).rx]',data(segment(i).lx :segment(i+1).rx),0);
      best = (0*( [segment(i).lx :segment(i+1).rx]' ))+coef(1);    
     segment(i).mc = sum((data([segment(i).lx :segment(i+1).rx]')-best).^2);
  end; 
          
end;


for i = 1 : length(segment) 
   
      coef = polyfit([segment(i).lx :segment(i).rx]',data(segment(i).lx :segment(i).rx),0);
      best = (0*( [segment(i).lx :segment(i).rx]' ))+coef(1);
      segment(i).y = best(1);
end;

residuals=[];
for i = 1 : length(segment)    
      coef = polyfit([segment(i).lx :segment(i).rx]',data(segment(i).lx :segment(i).rx),0);     
      best = (0*( [segment(i).lx :segment(i).rx]' ))+coef(1);
      residuals = [    residuals ; sum((data([segment(i).lx :segment(i).rx]')-best).^2)];
end;

if nargin == 3

hold on;
plot(data+0,'r');
print('apca','-dpng');
temp = [];

for i = 1 : length(segment) 
   
      coef = polyfit([segment(i).lx :segment(i).rx]',data(segment(i).lx :segment(i).rx),0);
      best = (0*( [segment(i).lx :segment(i).rx]' ))+coef(1);
      plot([segment(i).lx :segment(i).rx]', best,'b','LineWidth',3);
      temp = [temp; [best(1) best(end)]];    
end;

for i = 1 : length(segment)  - 1 
        plot([segment(i).rx :segment(i+1).lx]', [ temp(i,2) temp(i+1,1)  ],'g','LineWidth',3,'Color',[0 .7 0]);
end;
end


