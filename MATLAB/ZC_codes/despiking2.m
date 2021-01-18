function [spikes,newsignal,finalthreshold, window]=despiking2(signal) % injection peaks are replaced by a straight line joining the first and the last point of the peak

intsignal=signal;
newsignal=signal;

%sg2=[3,51,101,151,201,251,301,351,401,501,601,701,801,901,1001,1501];

dsignal=abs(diff(signal));

%Histogram of the diff signal
[counts,amp_range]=hist(dsignal,100);
counts=counts./length(dsignal);
threshold(1)=amp_range(2);

%Otsu thresholding method
%C1=find(amp_range<threshold(1));
proba_C1(1)=sum(counts(1));
proba_C2(1)=1-proba_C1(1);
mean_C1(1)=sum(amp_range(1).*counts(1))/proba_C1(1);
mean_C2(1)=sum(amp_range(2:length(amp_range)).*counts(2:length(amp_range)))/proba_C2(1);
for i=2:length(amp_range)-1
    threshold(i)=amp_range(i+1);
    %C1=find(amp_range<threshold(i));
    proba_C1(i)=sum(counts(1:i));
    proba_C2(i)=1-proba_C1(i);
    mean_C1(i)=sum(amp_range(1:i).*counts(1:i))/proba_C1(i);
    mean_C2(i)=sum(amp_range(i+1:length(amp_range)).*counts(i+1:length(amp_range)))/proba_C2(i);
    C1C2_var(i)=proba_C1(i)*proba_C2(i)*(mean_C1(i)-mean_C2(i))^2;
end
[max_C1C2var,k]=max(C1C2_var);
finalthreshold=threshold(k);

%iteration of the difference to find the window length
x=find((dsignal)> finalthreshold);
dx=diff(x);
y=find(dx>1);
dy=diff(y);
window=max(dy);
window=10*k;
window=window+1;

startspike=x(1);
i=1;;
for k=1:length(x)-1
      if x(k+1)-x(k)<=1
        k=k+1;
    else
        stopspike=x(k);
        spikes(i)=floor((stopspike+startspike)/2);
        startspike=x(k+1);
        
        % excision of the spike section in the original signal by linear
        % interpolation
        jmin(i)=max(1,spikes(i)-(window-1)/2); jmax(i)=min(spikes(i)+(window-1)/2,length(signal)); % since sgwindow is odd!
        l=[jmin(i):jmax(i)];
        intsignal(l)=intsignal(jmin(i))+(intsignal(jmax(i))-intsignal(jmin(i)))/((jmax(i)-jmin(i)))*(l-jmin(i)+1);
           
      
    i=i+1;
      end
    
end
         %last spike detection
        stopspike=x(length(x));
        spikes(i)=floor((stopspike+startspike)/2);
        jmin(i)=min(spikes(i)-(window-1)/2,length(signal)); jmax(i)=min(spikes(i)+(window-1)/2, length(signal));
        %excision of the last spike
        l=[jmin(i):jmax(i)];
        intsignal(l)=intsignal(jmin(i))+(intsignal(jmax(i))-intsignal(jmin(i)))/((jmax(i)-jmin(i)))*(l-jmin(i)+1);
         
        smooth=sgolayfilt(intsignal,2,3*window);
      
        for i=1:length(jmin)
            l=[jmin(i):jmax(i)];
            newsignal(l)=smooth(l);
        end
  % figure; plot(signal); hold on; plot(intsignal,'g'); plot(y1,'k'); plot(newsignal,'r');
  