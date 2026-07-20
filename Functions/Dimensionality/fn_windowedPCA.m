function WPCA = fn_windowedPCA(spike_conv, win_frames, ~)
% third argument accepted but unused — kept for API consistency
if nargin < 3, V_global = []; end   % third arg optional — ignored internally
[nN,nF]=size(spike_conv); nW=size(win_frames,1);
WPCA(nW)=struct(); CBLK=5000;
for w=1:nW
    f0=max(win_frames(w,1),1); f1=min(win_frames(w,2),nF);
    if f1<=f0, continue; end
    T_w=f1-f0+1; mu=mean(spike_conv(:,f0:f1),2); C=zeros(nN,nN);
    for b0=f0:CBLK:f1
        b1=min(b0+CBLK-1,f1); seg=double(spike_conv(:,b0:b1))-double(mu);
        C=C+seg*seg';
    end
    C=C/(T_w-1); [Vw,D]=eig(C); [ev,idx]=sort(diag(D),'descend');
    Vw=Vw(:,idx); ev(ev<0)=0;
    var_exp=100*ev/(sum(ev)+eps); cum_var=cumsum(var_exp);
    n95=find(cum_var>=95,1,'first'); if isempty(n95), n95=nN; end
    V3=Vw(:,1:3); sc=zeros(3,T_w,'single'); col=1;
    for b0=f0:CBLK:f1
        b1=min(b0+CBLK-1,f1); blk=b1-b0+1;
        seg=double(spike_conv(:,b0:b1))-double(mu);
        sc(:,col:col+blk-1)=single(V3'*seg); col=col+blk;
    end
    WPCA(w).win_frames=[f0,f1]; WPCA(w).t_axis=f0:f1;
    WPCA(w).eigenvectors=Vw; WPCA(w).eigenvalues=ev;
    WPCA(w).var_explained=var_exp; WPCA(w).cum_var=cum_var;
    WPCA(w).n_dims_95pct=n95; WPCA(w).scores3=sc;
    clear C Vw D seg sc V3 mu;
end
end
