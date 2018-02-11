function [R] = myinputEnter(St,D,Question)
% Input dialog. User may hit return after entering info.

if nargin<2 % Default values if called without args.
    D = '5';
elseif nargin<1
    St = 'myinput';
end

if ~ischar(D)
    D = num2str(D);
end

R = []; % In case the user closes the GUI.
S.fh = figure('units','pixels',...
              'position',[500 500 200 100],...
              'menubar','none',...
              'numbertitle','off',...
              'name',St,...
              'resize','off');
S.ed = uicontrol('style','edit',...
                 'units','pix',...
                'position',[10 60 180 30],...
                'string',D);
S.pb = uicontrol('style','pushbutton',...
                 'units','pix',...
                'position',[10 20 180 30],...
                'string',Question,...
                'callback',{@pb_call});
set(S.ed,'call',@ed_call)
uicontrol(S.ed) % Make the editbox active.
uiwait(S.fh) % Prevent all other processes from starting until closed.

    function [] = pb_call(varargin)
        R = str2double(get(S.ed,'string'));
        close(S.fh); % Closes the GUI, allows the new R to be returned.
    end

    function [] = ed_call(varargin)
        uicontrol(S.pb)
        drawnow
        R = str2double(get(S.ed,'string'));
        close(gcbf)
    end
end