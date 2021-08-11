function [firings]=ConvertAPtoST(input_firings,Pnrn)

function okin=MM_kin(x,y,n)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Converting Action potential to spike train
% input_firings = times of action action potential above certain threshold
% Pnrn = number of neurons

%% CODE

% Pnrn=16;
firings=[];
for i=1:Pnrn
    inds=input_firings((input_firings(:,2)==i));
    zero=zeros(1000,1);
    inds1=[zero;inds];
    isi=diff(inds1);
    ind=find(isi>10)-999;
    spikes=inds(ind);
    firings=[firings; spikes,i+0*spikes];
end

% [Y,I]=sort(firings(:,1));
% firings=firings(I,:);