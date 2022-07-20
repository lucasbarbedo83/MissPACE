function varargout = lisstvsf_view(varargin)
% LISSTVSF_VIEW MATLAB code for lisstvsf_view.fig
%      LISSTVSF_VIEW, by itself, creates a new LISSTVSF_VIEW or raises the existing
%      singleton*.
%
%      H = LISSTVSF_VIEW returns the handle to a new LISSTVSF_VIEW or the handle to
%      the existing singleton*.
%
%      LISSTVSF_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LISSTVSF_VIEW.M with the given input arguments.
%
%      LISSTVSF_VIEW('Property','Value',...) creates a new LISSTVSF_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lisstvsf_view_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lisstvsf_view_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lisstvsf_view

% Last Modified by GUIDE v2.5 20-Dec-2016 10:24:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lisstvsf_view_OpeningFcn, ...
                   'gui_OutputFcn',  @lisstvsf_view_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before lisstvsf_view is made visible.
function lisstvsf_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lisstvsf_view (see VARARGIN)

disp(varargin)


% Deal with issue where opengl rendering (with alphadata) mirrors tick
% labels
% opengl('software');
% % opengl VERBOSE
% % opengl INFO
% % if strfind(gl.Vendor,'ATI')
% %     opengl('software')
% % end

if ~exist('lisstvsf_settings.mat','file')
    Make_Settings;
end

handles.settings = load('lisstvsf_settings.mat');
handles.filename = '';
handles.comment = '';
handles.dat = [];
handles.comparedat = [];
set(handles.OverlayData,'enable','off')

% Plot Sequoia logo to figure 
set(handles.axes_logo, 'Visible','on');
axes(handles.axes_logo)
image(imread('sequoia_logo.jpg'));
axis off
handles.hzoom = zoom;
setAllowAxesZoom(handles.hzoom,handles.axes_logo,false);

set(handles.popup_plottype, 'String', {'Net Raw PMT Signals', 'Chopped Raw PMT Signals', ...
        'Ring Detector Scattering, LP, LREF', 'Auxiliary Parameters'});
set(handles.popup_plottype, 'Value', handles.settings.view_plottype)    

set(handles.uitoggletool_ToggleRaw, 'State',handles.settings.view_showraw);

% 
% if nargin>3
%     handles.mode = varargin{1};
%     handles.hMakep = varargin{2};
% else
%     handles.mode = 'normal';
%     handles.hMakep = [];
% end
% 
% if strcmpi(handles.mode, 'GetScat') || strcmpi(handles.mode, 'GetBackground')
%     
%     if nargin==6
%         appdata = getappdata(varargin{2},'AppData');
%     else
%         error('GetScat or GetBackground calling mode must have three arguments!')
%     end
%     
% %     set(hObject, 'WindowStyle','Modal')
%     
%     if strcmpi(handles.mode, 'GetScat')
%         set(hObject, 'Name','Select/Edit Scattering Data')
%     elseif strcmpi(handles.mode, 'GetBackground')
%         set(hObject, 'Name','Select/Edit Background Data')
%     end
%         
%     if isempty(varargin{3}) % dat structure was not provided in varargin
% 
%         [filename, pathname] = uigetfile({'*.DAT;*.VSF','All LISST-VSF Binary Data Files (*.DAT;*.VSF)'}, ...
%                                     'Open LISST-VSF Data File', ...
%                                     handles.settings.current_path);
%         if isequal(filename,0) || isequal(pathname,0)
%             disp('User pressed cancel: input file not selected')
%         else
%             handles.filename = fullfile(pathname, filename);
%             set(handles.text_inputfile, 'String',handles.filename); 
%             disp(['User selected: handles.filename = ' handles.filename])
%             handles.settings.current_path = pathname;
%             handles.dat = lisstvsf_readdat(handles.filename);
%             Refresh_Plot_Data(handles);
%         end
%     
%     else
%         
%         handles.dat = varargin{3};
%         set(handles.text_inputfile, 'String',['Previously loaded data from ' handles.dat.filename]); 
% 
%     end
% 
% end
    
Refresh_Plot_Data(handles)







% Choose default command line output for lisstvsf_view
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lisstvsf_view wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lisstvsf_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('figure1_DeleteFcn')

% Save application settings on close
disp('Saving application settings...');
settings = handles.settings;

[settings_path,tmp,tmp] = fileparts(mfilename('fullpath'));
%save([settings_path '\lisstvsf_settings.mat'],'-struct','settings');
save(fullfile(settings_path,'lisstvsf_settings.mat'),'-struct','settings');


% if ~isempty(handles.hMakep)
%     appdata = getappdata(handles.hMakep,'AppData');
%     if strcmpi(handles.mode, 'GetScat') 
%         appdata.dat_scat = handles.dat;
%     elseif strcmpi(handles.mode, 'GetBackground')
%         appdata.dat_zscat = handles.dat;
%     end
%     setappdata(handles.hMakep,'AppData', appdata);
% end
% -------------------------------------------------------------------------




% =========================================================================
% popup_plottype
% =========================================================================

% --- Executes on selection change in popup_plottype.
function popup_plottype_Callback(hObject, eventdata, handles)
% hObject    handle to popup_plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_plottype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_plottype
contents = cellstr(get(hObject,'String'));
get(hObject,'Value')
contents{get(hObject,'Value')};

handles.settings.view_plottype = get(hObject,'Value');

Refresh_Plot_Data(handles)


% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
% Helper to refresh data plots
% -------------------------------------------------------------------------
function Refresh_Plot_Data(handles)
disp('Setup_Plot_Axes')
plot_type = get(handles.popup_plottype, 'Value');

cla(handles.axes1)
cla(handles.axes2)
cla(handles.axes3)
cla(handles.axes4)

axis(handles.axes1, 'auto')
axis(handles.axes2, 'auto')
axis(handles.axes3, 'auto')
axis(handles.axes4, 'auto')

set(handles.axes1, 'Box','on', 'NextPlot','add')
set(handles.axes2, 'Box','on', 'NextPlot','add')
set(handles.axes3, 'Box','on', 'NextPlot','add')
set(handles.axes4, 'Box','on', 'NextPlot','add')

% Setup_Plot_Axes(handles);

ms = 5;
lw = 0.75;
lineColor = {'k','r','b','g'};

% {'Net PMT Signals', 'Ring Detector Scattering', ...
%          'Laser on/off PMT Signals', 'Aux Parameters'});
for n = 1:length(handles.dat)
switch plot_type

    % -----------------------------------------------------------------
    case 1, % Net PMT Signals
    % -----------------------------------------------------------------

        set(handles.axes4, 'Visible','on');

        set(get(handles.axes1,'XLabel'),'String','Angle [idx]')
        set(get(handles.axes2,'XLabel'),'String','Angle [idx]')
        set(get(handles.axes3,'XLabel'),'String','Angle [idx]')
        set(get(handles.axes4,'XLabel'),'String','Angle [idx]')
        set(get(handles.axes1,'YLabel'),'String','Net rp [counts]')
        set(get(handles.axes2,'YLabel'),'String','Net rr [counts]')
        set(get(handles.axes3,'YLabel'),'String','Net pp [counts]')
        set(get(handles.axes4,'YLabel'),'String','Net pr [counts]')
        set(handles.axes1, 'XScale','linear', 'YScale','log')
        set(handles.axes2, 'XScale','linear', 'YScale','log')
        set(handles.axes3, 'XScale','linear', 'YScale','log')
        set(handles.axes4, 'XScale','linear', 'YScale','log')

%         set([handles.axes1 handles.axes2 handles.axes3 handles.axes4], ...
%             'XLimMode','manual', 'XLim',[0 180])

        if ~isempty(handles.dat)
            if strcmp(handles.settings.view_showraw,'on')            
                plot(handles.axes1, handles.dat(n).angles_idx, handles.dat(n).rp, lineColor{n}, 'LineStyle', '--', 'LineWidth',lw);
                plot(handles.axes2, handles.dat(n).angles_idx, handles.dat(n).rr, lineColor{n}, 'LineStyle', '--', 'LineWidth',lw);
                plot(handles.axes3, handles.dat(n).angles_idx, handles.dat(n).pp, lineColor{n}, 'LineStyle', '--', 'LineWidth',lw);
                plot(handles.axes4, handles.dat(n).angles_idx, handles.dat(n).pr, lineColor{n}, 'LineStyle', '--', 'LineWidth',lw);
            end
            plot(handles.axes1, handles.dat(n).angles_idx, median(handles.dat(n).rp), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw*1.5);
            plot(handles.axes2, handles.dat(n).angles_idx, median(handles.dat(n).rr), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw*1.5);
            plot(handles.axes3, handles.dat(n).angles_idx, median(handles.dat(n).pp), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw*1.5);
            plot(handles.axes4, handles.dat(n).angles_idx, median(handles.dat(n).pr), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw*1.5);

        end

        legend(handles.axes1, 'off')
        legend(handles.axes2, 'off')
        legend(handles.axes3, 'off')
        legend(handles.axes4, 'off')
%         set([handles.axes1 handles.axes2 handles.axes3 handles.axes4], ...
%             'XLimMode','manual', 'XLim',[0 180])

%         zoom out
%         zoom reset

    % -----------------------------------------------------------------
    case 2, % Laser on/off PMT Signals
    % -----------------------------------------------------------------

        set(handles.axes4, 'Visible','on');

        set(get(handles.axes1,'XLabel'),'String','Angle [idx]')
        set(get(handles.axes2,'XLabel'),'String','Angle [idx]')
        set(get(handles.axes3,'XLabel'),'String','Angle [idx]')
        set(get(handles.axes4,'XLabel'),'String','Angle [idx]')
        set(get(handles.axes1,'YLabel'),'String','Laser on/off rp [counts]')
        set(get(handles.axes2,'YLabel'),'String','Laser on/off rr [counts]')
        set(get(handles.axes3,'YLabel'),'String','Laser on/off pp [counts]')
        set(get(handles.axes4,'YLabel'),'String','Laser on/off pr [counts]')
        set(handles.axes1, 'XScale','linear', 'YScale','linear')
        set(handles.axes2, 'XScale','linear', 'YScale','linear')
        set(handles.axes3, 'XScale','linear', 'YScale','linear')
        set(handles.axes4, 'XScale','linear', 'YScale','linear')

        if ~isempty(handles.dat(n))

           if strcmp(handles.settings.view_showraw,'on')            
                h1=plot(handles.axes1, handles.dat(n).angles_idx, handles.dat(n).raw.rp_off, lineColor{n}, 'LineStyle', '--', 'LineWidth',lw);
                h2=plot(handles.axes1, handles.dat(n).angles_idx, handles.dat(n).raw.rp_on, lineColor{n}, 'LineStyle', '-', 'LineWidth',lw);
                plot(handles.axes2, handles.dat(n).angles_idx, handles.dat(n).raw.rr_off, lineColor{n}, 'LineStyle', '--', 'LineWidth',lw);
                plot(handles.axes2, handles.dat(n).angles_idx, handles.dat(n).raw.rr_on, lineColor{n}, 'LineStyle', '-', 'LineWidth',lw);
                plot(handles.axes3, handles.dat(n).angles_idx, handles.dat(n).raw.pp_off, lineColor{n}, 'LineStyle', '--', 'LineWidth',lw);
                plot(handles.axes3, handles.dat(n).angles_idx, handles.dat(n).raw.pp_on, lineColor{n}, 'LineStyle', '-', 'LineWidth',lw);
                plot(handles.axes4, handles.dat(n).angles_idx, handles.dat(n).raw.pr_off, lineColor{n}, 'LineStyle', '--', 'LineWidth',lw);
                plot(handles.axes4, handles.dat(n).angles_idx, handles.dat(n).raw.pr_on, lineColor{n}, 'LineStyle', '-', 'LineWidth',lw);
            else
                h1=plot(handles.axes1, handles.dat(n).angles_idx, median(handles.dat(n).raw.rp_off), lineColor{n}, 'LineStyle', '--', 'LineWidth',lw*1.5);
                h2=plot(handles.axes1, handles.dat(n).angles_idx, median(handles.dat(n).raw.rp_on), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw*1.5);
                plot(handles.axes2, handles.dat(n).angles_idx, median(handles.dat(n).raw.rr_off), lineColor{n}, 'LineStyle', '--', 'LineWidth',lw*1.5);
                plot(handles.axes2, handles.dat(n).angles_idx, median(handles.dat(n).raw.rr_on), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw*1.5);
                plot(handles.axes3, handles.dat(n).angles_idx, median(handles.dat(n).raw.pp_off), lineColor{n}, 'LineStyle', '--', 'LineWidth',lw*1.5);
                plot(handles.axes3, handles.dat(n).angles_idx, median(handles.dat(n).raw.pp_on), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw*1.5);
                plot(handles.axes4, handles.dat(n).angles_idx, median(handles.dat(n).raw.pr_off), lineColor{n}, 'LineStyle', '--', 'LineWidth',lw*1.5);
                plot(handles.axes4, handles.dat(n).angles_idx, median(handles.dat(n).raw.pr_on), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw*1.5);
            end
            
            hl=legend(handles.axes1, [h1(1) h2(1)], {'Laser Off','Laser On'}, 'Location','best');
            set(hl,'Box','off', 'FontSize',round(get(handles.axes3,'FontSize')*0.8))
        end
        
        legend(handles.axes2, 'off')
        legend(handles.axes3, 'off')
        legend(handles.axes4, 'off')
        axis auto
%         set([handles.axes1 handles.axes2 handles.axes3 handles.axes4], ...
%             'XLimMode','manual', 'XLim',[0 180])
%         zoom out
%         zoom reset

    % -----------------------------------------------------------------
    case 3, % Ring Detector Scattering
    % -----------------------------------------------------------------

        set(handles.axes4, 'Visible','on');

        set(get(handles.axes1,'XLabel'),'String','Detector')
        set(get(handles.axes2,'XLabel'),'String','Detector')
        set(get(handles.axes3,'XLabel'),'String','Sample')
        set(get(handles.axes4,'XLabel'),'String','Sample')
        set(get(handles.axes1,'YLabel'),'String','Ring Scattering [counts]')
        set(get(handles.axes2,'YLabel'),'String','Ring Scattering [counts]')
        set(get(handles.axes3,'YLabel'),'String','Laser LP, LREF [counts]')
        set(get(handles.axes4,'YLabel'),'String','Laser LP/LREF Ratio')
        set(handles.axes1, 'XScale','linear', 'YScale','linear')
        set(handles.axes2, 'XScale','linear', 'YScale','linear')
        set(handles.axes3, 'XScale','linear', 'YScale','linear')
        set(handles.axes4, 'XScale','linear', 'YScale','linear')

        if ~isempty(handles.dat(n))
    
            if strcmp(handles.settings.view_showraw,'on')            
                plot(handles.axes1, 1:32, handles.dat(n).rings1, lineColor{n} ,'Marker', '.', 'MarkerSize',ms);
                plot(handles.axes2, 1:32, handles.dat(n).rings2, lineColor{n} ,'Marker', '.', 'MarkerSize',ms);
            end
            plot(handles.axes1, 1:32, median(handles.dat(n).rings1), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw);
            plot(handles.axes2, 1:32, median(handles.dat(n).rings2), lineColor{n}, 'LineStyle', '-', 'LineWidth',lw);

            h1=plot(handles.axes3, 1:size(handles.dat(n).LREF,1), handles.dat(n).LREF(:,1), lineColor{n},'Marker', 'o', 'LineStyle', '--', 'LineWidth',lw, 'MarkerSize',ms);
            h2=plot(handles.axes3, 1:size(handles.dat(n).LREF,1), handles.dat(n).LREF(:,2), lineColor{n},'Marker','s', 'LineStyle', '-', 'LineWidth',lw, 'MarkerSize',ms);
            h3=plot(handles.axes3, 1:size(handles.dat(n).LP,1), handles.dat(n).LP(:,1), lineColor{n},'Marker','d', 'LineStyle', '--', 'LineWidth',lw, 'MarkerSize',ms);
            h4=plot(handles.axes3, 1:size(handles.dat(n).LP,1), handles.dat(n).LP(:,2), lineColor{n},'Marker','^', 'LineStyle', '-', 'LineWidth',lw, 'MarkerSize',ms);
        
            hl=legend([h1(1) h2(1) h3(1) h4(1)], {'LREF (laser r)','LREF (laser p)','LP (laser r)','LP (laser p)'}, 'Location','best');
            set(hl,'Box','off', 'FontSize',round(get(handles.axes3,'FontSize')*0.8))
            
            h5=plot(handles.axes4, 1:size(handles.dat(n).LP,1), handles.dat(n).LP(:,1)./handles.dat(n).LREF(:,1), lineColor{n}, 'Marker' ,'o', 'LineStyle', '--', 'LineWidth',lw, 'MarkerSize',ms);
            h6=plot(handles.axes4, 1:size(handles.dat(n).LP,1), handles.dat(n).LP(:,2)./handles.dat(n).LREF(:,2), lineColor{n}, 'Marker', 'd', 'LineStyle', '-', 'LineWidth',lw, 'MarkerSize',ms);
 
            hl=legend([h5(1) h6(1)], {'Laser r','Laser p'}, 'Location','best');
            set(hl,'Box','off', 'FontSize',round(get(handles.axes3,'FontSize')*0.8))
        end

        legend(handles.axes1, 'off')
        legend(handles.axes2, 'off')
%         set([handles.axes1 handles.axes2],'XLimMode','manual', 'XLim',[0 33])
%         
%         set([handles.axes3 handles.axes4],'XLimMode','auto')
% %         zoom out
%         zoom reset

    % -----------------------------------------------------------------
    case 4, % Aux Parameters
    % -----------------------------------------------------------------

        set(get(handles.axes1,'XLabel'),'String','Sample')
        set(get(handles.axes2,'XLabel'),'String','Sample')
        set(get(handles.axes3,'XLabel'),'String','Sample')
        set(get(handles.axes4,'XLabel'),'String','Sample')
        set(get(handles.axes1,'YLabel'),'String','Depth [m]')
        set(get(handles.axes2,'YLabel'),'String','Temperature (slow) [^\circC]')
        set(get(handles.axes3,'YLabel'),'String','Battery [volts]')
        set(get(handles.axes4,'YLabel'),'String','PMT Gain')
        
        set(handles.axes1, 'XScale','linear', 'YScale','linear')
        set(handles.axes2, 'XScale','linear', 'YScale','linear')
        set(handles.axes3, 'XScale','linear', 'YScale','linear')
        set(handles.axes4, 'XScale','linear', 'YScale','linear')

        if ~isempty(handles.dat(n))
        
            %text(0.1,0.9, {'Filename', 'PMT Settings', 'Date and Time'}, ...
              %  'Parent',handles.axes4, ...
              %  'HorizontalAlignment','left', 'VerticalAlignment','top', ...
              %  'LineStyle','-');

            plot(handles.axes1, 1:length(handles.dat(n).depth), handles.dat(n).depth./100, lineColor{n}, 'Marker', '.', 'LineStyle', '-');
            plot(handles.axes2, 1:length(handles.dat(n).tempC), handles.dat(n).tempC./100, lineColor{n}, 'Marker', '.', 'LineStyle', '-');
            plot(handles.axes3, 1:length(handles.dat(n).batt_volts), handles.dat(n).batt_volts./100, lineColor{n}, 'Marker', '.', 'LineStyle', '-');
            plot(handles.axes4, 1:length(handles.dat(n).pmt_gain), handles.dat(n).pmt_gain, lineColor{n}, 'Marker', '.', 'LineStyle', '-');
        end
        legend(handles.axes1, 'off')
        legend(handles.axes2, 'off')
        legend(handles.axes3, 'off')
        legend(handles.axes4, 'off')
        
%         set([handles.axes1 handles.axes2 handles.axes3],'XLimMode','auto')
% %         zoom out
%         zoom reset


end
end








% =========================================================================
% edit_comment
% =========================================================================
function edit_comment_Callback(hObject, eventdata, handles)
% hObject    handle to edit_comment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_comment as text
%        str2double(get(hObject,'String')) returns contents of edit_comment as a double

handles.comment = get(hObject,'String');


% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_comment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_comment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% =========================================================================
% Pushbutton copy
% =========================================================================
function uipushtool_copy_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hnewfig = figure('Visible','off');
hax1 = copyobj(handles.axes1, hnewfig);
hax2 = copyobj(handles.axes2, hnewfig);
hax3 = copyobj(handles.axes3, hnewfig);
hax4 = copyobj(handles.axes4, hnewfig);

left1 = 0.1;
left2 = 0.6;
bottom1 = 0.6;
bottom2 = 0.1;
width = 0.38;
height = 0.38;

%[left bottom width height]
set(hax1, 'Units','normalized', 'Position',[left1 bottom1 width height])
set(hax2, 'Units','normalized', 'Position',[left2 bottom1 width height])
set(hax3, 'Units','normalized', 'Position',[left1 bottom2 width height])
set(hax4, 'Units','normalized', 'Position',[left2 bottom2 width height])

% set(0,'showhiddenhandles','on')
% set(hnewfig, 'PaperPositionMode','auto')
print(hnewfig, '-dmeta')


% =========================================================================
% Pushbutton FileOpen
% =========================================================================
function uipushtool_FileOpen_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_FileOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %filespec = fullfile(handles.settings.current_path)
    %     handles.settings

    [filename, pathname] = uigetfile({'*.DAT;*.VSF','All LISST-VSF Binary Data Files (*.DAT;*.VSF)'}, ...
                                'Open LISST-VSF Data File', ...
                                handles.settings.current_path);
                            
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel: input file not selected')
    else
        handles.filename = fullfile(pathname, filename);
        set(handles.text_inputfile, 'String',handles.filename); 
        disp(['User selected: handles.filename = ' handles.filename])
        handles.settings.current_path = pathname;
        handles.settings
        
        handles.dat = lisstvsf_readdat(handles.filename);
        handles.dat
        
        Refresh_Plot_Data(handles);
        
        
        
    end
    
% Update handles structure
guidata(hObject, handles);


% =========================================================================
% Pushbutton ExportToWorkspace
% =========================================================================
function uipushtool_ExportToWorkspace_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_ExportToWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.dat)
    
    handles.dat(1).comment = handles.comment;
    
    [tmp,fn] = fileparts(handles.filename);
    fn = strrep(fn, '#', '_');
    fn = strtok(fn, '.');
    assignin('base', fn, handles.dat);
end


% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function uitoggletool_ToggleRaw_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_ToggleRaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.settings.view_showraw = get(hObject, 'State');

Refresh_Plot_Data(handles)

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function uipushtool_Save_ClickedCallback(hObject, eventdata, handles)

    [filename, pathname] = uiputfile({'*.png','png files (*.png)'}, ...
                                'Save window as a png file', ...
                                handles.settings.current_path);
 print(fullfile(pathname,filename),'-dpng')


% --- Executes on button press in pushbutton_FileOpen.
function pushbutton_FileOpen_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_FileOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [filename, pathname] = uigetfile({'*.DAT;*.VSF','All LISST-VSF Binary Data Files (*.DAT;*.VSF)'}, ...
                                'Open LISST-VSF Data File', ...
                                handles.settings.current_path);
                            
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel: input file not selected')
    else
        handles.filename = fullfile(pathname, filename);
        set(handles.text_inputfile, 'String',handles.filename); 
        disp(['User selected: handles.filename = ' handles.filename])
        handles.settings.current_path = pathname;
        handles.settings
        
        handles.dat = lisstvsf_readdat(handles.filename);
        handles.dat
        
        Refresh_Plot_Data(handles);
        
        set(handles.OverlayData,'enable','on')
        handles.allFilenames = handles.filename;
    end
    
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_showfzsc.
function checkbox_showfzsc_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showfzsc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showfzsc


% --- Executes on button press in OverlayData.
function OverlayData_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile({'*.DAT;*.VSF','All LISST-VSF Binary Data Files (*.DAT;*.VSF)'}, ...
                                'Open LISST-VSF Data File', ...
                                handles.settings.current_path);
                            
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel: input file not selected')
    else
        handles.allFilenames = [handles.allFilenames char(10) fullfile(pathname, filename)]; 
        handles.filename = fullfile(pathname, filename);
        set(handles.text_inputfile, 'String',handles.allFilenames); 
        disp(['User selected: handles.filename = ' handles.allFilenames])
        handles.settings.current_path = pathname;
        handles.settings
        
        handles.dat(end+1) = lisstvsf_readdat(handles.filename);
        
        Refresh_Plot_Data(handles);
        
        if length(handles.dat) > 3
            set(handles.OverlayData,'enable','off')
        end
        
    end
    
% Update handles structure
guidata(hObject, handles);

function Make_Settings

current_path = '';

view_plottype = 1;
view_showraw = 'on';

makep_fzsfile = '';
makep_calfactors = '';
makep_plottype = 1;
makep_showraw = 'on';

filename_fzsc = '';
filename_calfact = '';

save lisstvsf_settings.mat
