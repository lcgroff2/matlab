
% sifreport - report on kinetics traces from sifkin

% first show particle positions

figure(2)
hold off
imagesc(datacrop.image(:,:,1)'), shading flat
set(gca,'dataaspectratio',[1 1 1])
colorbar
hold on % keep image onscreen while markers are added

for idx = 1:length(pbkin.xpos);
  xpos=pbkin.xpos{idx};
  ypos=pbkin.ypos{idx};
  plot(xpos,ypos,'go')
  txth=text(xpos+n+2,ypos,num2str(idx));
  set(txth,'color',[0 1 0]); % color green
  rectangle('position',[xpos-n,ypos-n,2*n,2*n],'edgecolor','r')
  fprintf(1,'Particle %i: Deathnum = %1.4g, Tau_est = %1.4g s, Decay %2.3g pct\n',idx,pbkin.deathnum{idx},pbkin.tau{idx},pbkin.frac_decay{idx}*100);
end
hold off % next plot gets a blank slate
