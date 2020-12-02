close all; clc; clear all;

fig1 = openfig('EF1.fig');
fig2 = openfig('EF2.fig');
fig3 = openfig('EF3.fig');
fig4 = openfig('EF4.fig');
fig5 = openfig('EF5.fig');
fig6 = openfig('EF6.fig');

fig2.Children.Position = fig1.Children.Position;
fig3.Children.Position = fig1.Children.Position;
fig4.Children.Position = fig1.Children.Position;
fig5.Children.Position = fig1.Children.Position;
fig6.Children.Position = fig1.Children.Position;

fig2.Children.XTickLabel = fig1.Children.XTickLabel;
fig3.Children.XTickLabel = fig1.Children.XTickLabel;
fig4.Children.XTickLabel = fig1.Children.XTickLabel;
fig5.Children.XTickLabel = fig1.Children.XTickLabel;
fig6.Children.XTickLabel = fig1.Children.XTickLabel;

fig2.Children.TickLabelInterpreter = fig1.Children.TickLabelInterpreter;
fig3.Children.TickLabelInterpreter = fig1.Children.TickLabelInterpreter;
fig4.Children.TickLabelInterpreter = fig1.Children.TickLabelInterpreter;
fig5.Children.TickLabelInterpreter = fig1.Children.TickLabelInterpreter;
fig6.Children.TickLabelInterpreter = fig1.Children.TickLabelInterpreter;

fig2.Children.YGrid = 'on';
fig3.Children.YGrid = 'on';
fig4.Children.YGrid = 'on';
fig5.Children.YGrid = 'on';
fig6.Children.YGrid = 'on';

fig2.GraphicsSmoothing = 'off';
fig3.GraphicsSmoothing = 'off';
fig4.GraphicsSmoothing = 'off';
fig5.GraphicsSmoothing = 'off';
fig6.GraphicsSmoothing = 'off';

saveas(fig1,'EF1_boxplot_media_mediana.png')
saveas(fig2,'EF2_boxplot_media_mediana.png')
saveas(fig3,'EF3_boxplot_media_mediana.png')
saveas(fig4,'EF4_boxplot_media_mediana.png')
saveas(fig5,'EF5_boxplot_media_mediana.png')
saveas(fig6,'EF6_boxplot_media_mediana.png')