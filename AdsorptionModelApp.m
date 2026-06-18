classdef AdsorptionModelApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure               matlab.ui.Figure
        GridLayout             matlab.ui.container.GridLayout
        LeftPanel              matlab.ui.container.Panel
        RunButton              matlab.ui.control.Button
        ParameterPanel         matlab.ui.container.Panel
        ColumnLengthEditField  matlab.ui.control.NumericEditField
        ColumnLengthLabel      matlab.ui.control.Label
        SpatialPointsEditField matlab.ui.control.NumericEditField
        SpatialPointsLabel     matlab.ui.control.Label
        VelocityEditField      matlab.ui.control.NumericEditField
        VelocityLabel          matlab.ui.control.Label
        DispersionEditField    matlab.ui.control.NumericEditField
        DispersionLabel        matlab.ui.control.Label
        PorosityEditField      matlab.ui.control.NumericEditField
        PorosityLabel          matlab.ui.control.Label
        DensityEditField       matlab.ui.control.NumericEditField
        DensityLabel           matlab.ui.control.Label
        qmaxEditField          matlab.ui.control.NumericEditField
        qmaxLabel              matlab.ui.control.Label
        kaEditField            matlab.ui.control.NumericEditField
        kaLabel                matlab.ui.control.Label
        kdEditField            matlab.ui.control.NumericEditField
        kdLabel                matlab.ui.control.Label
        C0EditField            matlab.ui.control.NumericEditField
        C0Label                matlab.ui.control.Label
        TimeSpanEditField      matlab.ui.control.NumericEditField
        TimeSpanLabel          matlab.ui.control.Label
        RightPanel             matlab.ui.container.Panel
        ResultsTabGroup        matlab.ui.container.TabGroup
        ResultsTab             matlab.ui.container.Tab
        PlotAxes1              matlab.ui.control.UIAxes % Top-left
        PlotAxes2              matlab.ui.control.UIAxes % Top-right
        PlotAxes3              matlab.ui.control.UIAxes % Bottom-left
        PlotAxes4              matlab.ui.control.UIAxes % Bottom-right
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: RunButton
        function RunButtonPushed(app, ~)
            try
                % Get all parameters from UI
                L = app.ColumnLengthEditField.Value;
                n = app.SpatialPointsEditField.Value;
                u = app.VelocityEditField.Value;
                D_ax = app.DispersionEditField.Value;
                eps = app.PorosityEditField.Value;
                rho_s = app.DensityEditField.Value;
                q_max = app.qmaxEditField.Value;
                k_a = app.kaEditField.Value;
                k_d = app.kdEditField.Value;
                C_in = app.C0EditField.Value;
                t_final = app.TimeSpanEditField.Value;

                % Validate inputs
                if L <= 0 || n < 2 || u <= 0 || D_ax < 0 || eps <= 0 || eps >= 1 || ...
                   rho_s <= 0 || q_max <= 0 || k_a < 0 || k_d < 0 || C_in < 0 || t_final <= 0
                    uialert(app.UIFigure, 'Invalid parameters: Ensure positive values and 0 < ε < 1.', 'Input Error');
                    return;
                end

                % Calculate spatial step
                dz = L / (n-1);

                % Initial conditions
                C0 = zeros(n,1);
                q0 = zeros(n,1);
                y0 = [C0; q0];
                tspan = [0 t_final];

                % Solve ODE with updated function
                options = odeset('RelTol',1e-3,'AbsTol',1e-5);
                [t, y] = ode15s(@(t,y) packedBedODE(t,y,n,dz,u,D_ax,eps,rho_s,k_a,k_d,q_max,C_in), tspan, y0, options);

                % Extract solutions
                C = y(:,1:n);
                q = y(:,n+1:end);
                z = linspace(0, L, n);
                
                % Clear previous plots
                cla(app.PlotAxes1, 'reset');
                cla(app.PlotAxes2, 'reset');
                cla(app.PlotAxes3, 'reset');
                cla(app.PlotAxes4, 'reset');

                % Plot results in 2x2 layout
                plot(app.PlotAxes1, z, C(end,:), 'b-', 'LineWidth', 2);
                xlabel(app.PlotAxes1, 'Column length z (m)');
                ylabel(app.PlotAxes1, 'C(z,t_{end}) [mol/m³]');
                title(app.PlotAxes1, 'CO₂ Concentration Profile at Final Time');
                grid(app.PlotAxes1, 'on');

                semilogx(app.PlotAxes2, t, C(:,end), 'r-', 'LineWidth', 2);
                xlabel(app.PlotAxes2, 'Time (s)');
                ylabel(app.PlotAxes2, 'C(L,t) [mol/m³]');
                title(app.PlotAxes2, 'Breakthrough Curve at Column Outlet');
                grid(app.PlotAxes2, 'on');
                app.PlotAxes2.YAxis.TickLabelFormat = '%.1e';

                plot(app.PlotAxes3, z, q(end,:), 'g-', 'LineWidth', 2);
                xlabel(app.PlotAxes3, 'Column length z (m)');
                ylabel(app.PlotAxes3, 'q(z,t_{end}) [mol/kg]');
                title(app.PlotAxes3, 'Adsorbed Phase Profile at Final Time');
                grid(app.PlotAxes3, 'on');

                semilogx(app.PlotAxes4, t, q(:,end), 'm-', 'LineWidth', 2);
                xlabel(app.PlotAxes4, 'Time (s)');
                ylabel(app.PlotAxes4, 'q(L,t) [mol/kg]');
                title(app.PlotAxes4, 'Adsorption at Column Outlet Over Time');
                grid(app.PlotAxes4, 'on');
                app.PlotAxes4.YAxis.TickLabelFormat = '%.1e';

                % Ensure axes are visible
                app.PlotAxes1.Visible = 'on';
                app.PlotAxes2.Visible = 'on';
                app.PlotAxes3.Visible = 'on';
                app.PlotAxes4.Visible = 'on';
                drawnow; % Force plot update

            catch ME
                % Display error message
                uialert(app.UIFigure, ['Error: ' ME.message], 'Simulation Error');
                disp(ME.message); % Debug in command window
            end
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            % Create UIFigure with increased height
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 1200 800]; % Height remains 800
            app.UIFigure.Name = 'CO₂ Adsorption in Packed Bed (Optimized)';
            
            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x', '2x'};
            app.GridLayout.RowHeight = {'1x'};
            
            % Left Panel with Parameters
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;
            
            % Parameter Panel with adjusted position
            app.ParameterPanel = uipanel(app.LeftPanel);
            app.ParameterPanel.Title = 'Model Parameters (Optimized)';
            app.ParameterPanel.Position = [50 70 280 700]; % y-position at 70
            
            % Create parameter fields with default values
            createParameterField(app, 'ColumnLength', 'Column length (m)', 1.0, [20 620 240 22]); % Start at 640
            createParameterField(app, 'SpatialPoints', 'Spatial points (n)', 200, [20 563 240 22]); % 640 - 57 
            createParameterField(app, 'Velocity', 'Velocity u (m/s)', 0.01, [20 506 240 22]);     % 583 - 57
            createParameterField(app, 'Dispersion', 'Dispersion Dₐₓ (m²/s)', 1e-6, [20 449 240 22]); % 526 - 57
            createParameterField(app, 'Porosity', 'Porosity ε', 0.35, [20 392 240 22]);         % 469 - 57
            createParameterField(app, 'Density', 'Density ρₛ (kg/m³)', 50, [20 335 240 22]);    % 412 - 57
            createParameterField(app, 'qmax', 'Capacity qₘₐₓ (mol/kg)', 1.0, [20 278 240 22]);  % 355 - 57
            createParameterField(app, 'ka', 'Adsorption kₐ (1/s)', 0.01, [20 221 240 22]);      % 298 - 57
            createParameterField(app, 'kd', 'Desorption k_d (1/s)', 0.005, [20 164 240 22]);    % 241 - 57
            createParameterField(app, 'C0', 'Inlet C₀ (mol/m³)', 1.0, [20 107 240 22]);          % 184 - 57
            createParameterField(app, 'TimeSpan', 'Time span (s)', 5000, [20 50 240 22]);       % 127 - 57
            
            % Run Button (adjusted position for symmetry)
            app.RunButton = uibutton(app.LeftPanel, 'push');
            app.RunButton.Position = [130 30 120 30]; % Moved up to 80
            app.RunButton.Text = 'Run Simulation';
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            
            
            
            % Right Panel with Results
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            
            % Tab Group for Results
            app.ResultsTabGroup = uitabgroup(app.RightPanel);
            app.ResultsTabGroup.Position = [0 10 850 770];
            
            % Create single tab with 2x2 UIAxes
            app.ResultsTab = uitab(app.ResultsTabGroup, 'Title', 'Results');
            app.PlotAxes1 = uiaxes(app.ResultsTab);
            app.PlotAxes1.Position = [50 400 350 300]; % Top-left
            app.PlotAxes1.FontSize = 12;

            app.PlotAxes2 = uiaxes(app.ResultsTab);
            app.PlotAxes2.Position = [450 400 350 300]; % Top-right
            app.PlotAxes2.FontSize = 12;

            app.PlotAxes3 = uiaxes(app.ResultsTab);
            app.PlotAxes3.Position = [50 50 350 300]; % Bottom-left
            app.PlotAxes3.FontSize = 12;

            app.PlotAxes4 = uiaxes(app.ResultsTab);
            app.PlotAxes4.Position = [450 50 350 300]; % Bottom-right
            app.PlotAxes4.FontSize = 12;
        end
        
        function createParameterField(app, fieldName, labelText, defaultValue, position)
            % Create label with adjusted height and spacing
            label = uilabel(app.ParameterPanel);
            label.Text = labelText;
            label.Position = [position(1) position(2) + 35 position(3) 22]; % 35-pixel spacing, 22-pixel height
            label.FontWeight = 'bold';
            label.VerticalAlignment = 'bottom'; % Align text to bottom for consistency
            
            % Create edit field
            editField = uieditfield(app.ParameterPanel, 'numeric');
            editField.Value = defaultValue;
            editField.Position = position; % [x y width height]
            editField.RoundFractionalValues = 'off'; % Ensure decimal inputs
            
            app.([fieldName 'Label']) = label;
            app.([fieldName 'EditField']) = editField;
        end
    end

    methods (Access = public)
        % Construct app
        function app = AdsorptionModelApp
            createComponents(app);
            registerApp(app, app.UIFigure);
            if nargout == 0
                clear app;
            end
        end

        function delete(app)
            delete(app.UIFigure);
        end
    end
end

% Updated ODE function with improved numerical stability
function dydt = packedBedODE(~, y, n, dz, u, D_ax, eps, rho_s, k_a, k_d, q_max, C_in)
    % Split variables with physical bounds
    C = max(0, y(1:n));
    q = max(0, min(y(n+1:end), q_max*0.999)); % Prevent q > q_max
    
    % Initialize derivatives
    dCdt = zeros(n,1);
    dqdt = zeros(n,1);
    
    % Boundary conditions
    C_left = C_in;
    C_right = C(end); % Open boundary
    
    for i = 1:n
        % Upwind differencing for stability
        if i == 1
            dCdz = (C(1) - C_left)/dz;
            d2Cdz2 = (C(2) - 2*C(1) + C_left)/dz^2;
        elseif i == n
            dCdz = (C(n) - C(n-1))/dz;
            d2Cdz2 = (C_right - 2*C(n) + C(n-1))/dz^2;
        else
            dCdz = (C(i) - C(i-1))/dz;
            d2Cdz2 = (C(i+1) - 2*C(i) + C(i-1))/dz^2;
        end
        
        % Rate-limited adsorption kinetics
        q_avail = max(0, q_max - q(i));
        dqdt(i) = min(0.1, max(-0.1, k_a*C(i)*q_avail - k_d*q(i)));
        
        % Mass balance with rate limiting
        dCdt(i) = (1/eps)*(D_ax*d2Cdz2 - u*dCdz) - (1-eps)/eps*rho_s*dqdt(i);
    end
    
    dydt = [dCdt; dqdt];
end