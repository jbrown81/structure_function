function [] =  plotcorr(x,y,varargin)

% plotcorr.m
% Copyright (c) 2025 Radiata, Inc.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Given two input vectors:
% 1) regress y on x and return the F stat and p-value
% 2) sort the points in the x vector and plot x vs y
% 3) plot the regression line

plotcolor='blue';
point_labels={};
if nargin == 3, plotcolor = varargin{1}; end
if nargin == 4, plotcolor = varargin{1}; point_labels=varargin{2}; end

% transpose row vectors to column vectors
[xr,xc]=size(x);
[yr,yc]=size(y);
if xr==1
    x=x';
end
if yr==1
    y=y';
end

% regress y on x
[b,bint,r,rint,stats]=regress(y,[ones(length(x),1) x]);
f=stats(2);
p=stats(3);
%[b,stats]=robustfit(x,y);
%regp=stats.p;

[r,p]=corrcoef(y,x);
r=r(1,2);
p=p(1,2);
disp(strcat('correlation r-value: ',num2str(r)))
disp(strcat('correlation p-value: ',num2str(p)))

%plot points based on x sort
points=scatter(x,y,'filled',plotcolor);
set(points,'LineWidth',1,'Marker', 'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',plotcolor);
hold on;
axis tight
if ~isempty(point_labels)
    text(x(i1),y(i1),point_labels(i1));
end

% plot regression line
x_span=max(x)-min(x);
x_buffer=x_span*.1;
xreg=linspace(min(x)-x_buffer,max(x)+x_buffer,10);
yreg=b(1)+(b(2).*xreg);
regline=plot(xreg , yreg,'-');
set(regline,'LineWidth',2,'Color',plotcolor);
axis tight
legend(sprintf('r=%0.3f p=%0.3f',r,p),'Location','Best')

% plot with dif colors for group 1 vs. 2
%for i=1:40
%if i1(i)>20
%scatter(val1(i),ysort(i),'filled','red')
%else
%scatter(val1(i),ysort(i),'filled','blue')
%end
%end