% MATLAB script to batch-process AFM approach force curves stored in an .xlsx file
% Force curve data format (.xlsx):
% - Approach portion only, baseline flattened, calibrated
% - Excel layout: columns alternate: col1 = distance (nm), col2 = force (nN), col3 = distance (nm), col4 = force (nN), etc.
%
% Model: Bottom-effect (finite-thickness) correction applied using S. Chiodini et al., Bottom Effect in Atomic Force Microscopy Nanomechanics. Small 16,  (2020).
%   psi(chi) = 1 + 1.133*chi + 1.283*chi^2 + 0.769*chi^3 + 0.0975*chi^4
%   where chi = a/h and a = sqrt(R*delta) is the contact radius.
%
% Interactive workflow (per curve):
% 1) Show one force-vs-distance curve (force in nN on y, distance in nm on
%    x) at a time
% 2) User drags two vertical lines to define the fitting region.
% 3) Script fits E (Young's modulus) using corrected Hertz model on selected region
%    and overlays fitted curve and reports E on the plot.
% 4) User can choose to Save this curve and move to next curve or reselect/exit.

clear; close all; clc;

[filename, path] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx,*.xls)'}, 'Select AFM Excel file');
if isequal(filename,0), error('No file selected.'); end
fullfile_in = fullfile(path, filename);

prompt = {'Tip radius R (nm)','Poisson''s ratio \nu (use 0.5 for incompressible)','Film thickness h (nm)'};
defaults = {'5','0.5','4.3'};
answer = inputdlg(prompt,'Model parameters',1,defaults);
if isempty(answer), error('Operation cancelled by user.'); end
R_nm = str2double(answer{1});
nu = str2double(answer{2});
h_nm = str2double(answer{3});

R = R_nm*1e-9; % m
h = h_nm*1e-9; % m

[~,name_noext,~] = fileparts(filename);
out_allcurves_csv = fullfile(path, [name_noext '_raw_and_fits.csv']);
out_moduli_csv = fullfile(path, [name_noext '_moduli().csv']);

[data,~,~] = xlsread(fullfile_in);
if isempty(data), error('No numeric data found.'); end

[nrows, ncols] = size(data);
npairs = floor(ncols/2);

fixedSpacing_nm = 1;  % fixed selection range, unit: nm

fid = fopen(out_moduli_csv,'w');
fprintf(fid, 'curve_id, pair_index, E_Pa, E_MPa, R_nm, nu, h_nm, n_points_used, fit_R2\n');
fclose(fid);

curves = struct();
counter = 0;
for p = 1:npairs
    D = data(:,2*p-1);
    F = data(:,2*p);
    valid = ~isnan(F) & ~isnan(D);
    F = F(valid);
    D = D(valid);
    if numel(F)<10, continue; end
    counter = counter + 1;
    curves(counter).F_nN = F;
    curves(counter).D_nm = D;
    curves(counter).pair = p;
end

ncurves = numel(curves);
fprintf('Found %d curve(s) to process.\n', ncurves);
if ncurves==0, error('No valid curves found.'); end

idx = 1;
while idx < (ncurves + 1)
    F_nN = curves(idx).F_nN;
    D_nm = curves(idx).D_nm;
    F = F_nN*1e-9; % N
    D = D_nm*1e-9; % m

    hfig = figure('Name', sprintf('Curve %d/%d', idx, ncurves));
    plot(D*1e9,F*1e9,'-b'); hold on; grid on;
    xlabel('Distance (nm)'); ylabel('Force (nN)');
    title('Drag vertical lines to select fit region, then press Enter');

    % Create two draggable vertical lines
    yl = ylim;
    hLine1 = line([mean(D_nm)+7 mean(D_nm)+7], yl, 'Color','r','LineWidth',1.5,'ButtonDownFcn',@(src,~) dragLine(src,[],fixedSpacing_nm,hLine2));
    hLine2 = line([mean(D_nm)+7+fixedSpacing_nm mean(D_nm)+7+fixedSpacing_nm], yl, 'Color','r','LineWidth',1.5,'ButtonDownFcn',@(src,~) dragLine(src,[],fixedSpacing_nm,hLine1));


    % Wait for user to press Enter
    pause; % wait for key press

    x1 = get(hLine1,'XData'); x1 = x1(1);
    x2 = get(hLine2,'XData'); x2 = x2(1);
    if x2<x1, tmp=x1; x1=x2; x2=tmp; end

    sel_idx = find(D_nm>=x1 & D_nm<=x2);
    if numel(sel_idx)<5
        warning('Too few points selected, skipping curve.');
        close(hfig); continue;
    end

    delta = D(sel_idx) - D(sel_idx(1));
    Fsel = F(sel_idx) - F(sel_idx(1));

    sqrtR = sqrt(R);
    a = sqrt(R*delta);
    chi = a./h;
    psi = 1 + 1.133.*chi + 1.283.*chi.^2 + 0.769.*chi.^3 + 0.0975.*chi.^4;
    K = (4/3)/(1-nu^2)*sqrtR.*delta.^(3/2).*psi;

    E_est = (K'*Fsel)/(K'*K);
    F_fit = K*E_est + F(sel_idx(1));
    ss_res = sum((Fsel-F_fit+F(sel_idx(1))).^2);
    ss_tot = sum((Fsel-mean(Fsel)).^2);
    R2 = 1-ss_res/ss_tot;

    E_Pa = E_est; E_MPa = E_Pa/1e6;

    plot(D(sel_idx)*1e9,F_fit*1e9,'--r','LineWidth',1.5);
    legend('Raw','Fit');
    txt = sprintf('E = %.3g MPa (R^2=%.3f)',E_MPa,R2);
    text(0.05*max(D_nm),0.9*max(F_nN),txt,'FontSize',12,'BackgroundColor','w');

    % fid = fopen(out_allcurves_csv,'a');
    % fprintf(fid,'\n# Curve %d — pair %d\n',idx,curves(idx).pair);
    % fprintf(fid,'distance_nm, raw_force_nN, fit_force_nN\n');
    % for ii=1:numel(D)
    %     fprintf(fid,'%.6f, %.6f\n',D(ii)*1e9, F(ii)*1e9);
    % end
    % fclose(fid);

    choice = questdlg('Choose action','Next step','Next curve','Reselect','Exit','Next curve');
    if strcmp(choice,'Reselect')
        close(hfig); continue;
    elseif strcmp(choice,'Exit')
        close(hfig); break;
    else
        close(hfig);
        fid = fopen(out_moduli_csv,'a');
        fprintf(fid,'%d, %d, %.6e, %.6f, %.2f, %.3f, %.2f, %d, %.6f\n',...
            idx,curves(idx).pair,E_Pa,E_MPa,R_nm,nu,h_nm,numel(sel_idx),R2);
        fclose(fid);
        idx = idx + 1; 
        continue;
    end
end

fprintf('Processing complete. Results saved to:\n %s\n %s\n',out_allcurves_csv,out_moduli_csv);

function dragLine(src,~,fixedSpacing_nm,otherLine)
    fig = ancestor(src,'figure');
    set(fig,'WindowButtonMotionFcn',@(fig,~) moveLine(src,otherLine,fixedSpacing_nm));
    set(fig,'WindowButtonUpFcn',@(fig,~) stopDrag(fig));
end

function moveLine(lineObj,otherLine,fixedSpacing_nm)
    currPt = get(gca,'CurrentPoint');
    x = currPt(1,1);
    set(lineObj,'XData',[x x]);

    % judge which line is being dragged
    x_other = get(otherLine,'XData'); 
    x_other = x_other(1);

    % select a fixed range based on the dragged line
    if x < x_other
        newX = x + fixedSpacing_nm;
        set(otherLine,'XData',[newX newX]);
    else
        newX = x - fixedSpacing_nm;
        set(otherLine,'XData',[newX newX]);
    end
end


function stopDrag(fig)
set(fig,'WindowButtonMotionFcn','');
set(fig,'WindowButtonUpFcn','');
end
