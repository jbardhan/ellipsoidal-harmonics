figure;
set(gca,'fontsize',16);
kcalConv=332.112;

CFA_H = plot(Eexact,Eexact_cfa,'rx');
set(CFA_H,'MarkerSize',6,'linewidth',2);
hold on;
P_H = plot(Eexact,Eexact_p,'bo');
set(P_H,'MarkerSize',8,'linewidth',2);
C_H = plot(Eexact,Eexact_c,'ks');
set(C_H,'MarkerSize',10,'linewidth',2);
BEMC_H = plot(Eexact,Ebem_c,'md');
set(BEMC_H,'MarkerSize',4,'linewidth',1.5);
XY_H=plot([0 min(Eexact)],[0 min(Eexact)],'k');
set(XY_H,'Linewidth',2);

%axis([min(Eexact)*1.1 0 min(Eexact_p)*1.05 0]);
axis([-80 0 -140 0]);
ylabel('Estimated Solvation Free Energy (kcal/mol)');
xlabel('Analytical Solvation Free Energy (kcal/mol)');
LegH=legend('Exact BIBEE/CFA','Exact BIBEE/P','Exact BIBEE/M', 'BEM BIBEE/M','location','southeast');
print -depsc2 figure2-correction-impact.eps
print -dill figure2-correction-impact.ai
print -dpng figure2-correction-impact.png
