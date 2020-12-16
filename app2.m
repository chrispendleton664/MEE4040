classdef app2 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        GridLayout                  matlab.ui.container.GridLayout
        LeftPanel                   matlab.ui.container.Panel
        E1EditFieldLabel            matlab.ui.control.Label
        E1EditField                 matlab.ui.control.NumericEditField
        E2EditFieldLabel            matlab.ui.control.Label
        E2EditField                 matlab.ui.control.NumericEditField
        v12EditFieldLabel           matlab.ui.control.Label
        v12EditField                matlab.ui.control.NumericEditField
        G12EditFieldLabel           matlab.ui.control.Label
        G12EditField                matlab.ui.control.NumericEditField
        fxEditFieldLabel            matlab.ui.control.Label
        fxEditField                 matlab.ui.control.NumericEditField
        fyEditFieldLabel            matlab.ui.control.Label
        fyEditField                 matlab.ui.control.NumericEditField
        fxyEditFieldLabel           matlab.ui.control.Label
        fxyEditField                matlab.ui.control.NumericEditField
        thetaEditFieldLabel         matlab.ui.control.Label
        thetaEditField              matlab.ui.control.NumericEditField
        thetaEditField_2Label       matlab.ui.control.Label
        thetaEditField_2            matlab.ui.control.NumericEditField
        thetaEditField_3Label       matlab.ui.control.Label
        thetaEditField_3            matlab.ui.control.NumericEditField
        NoPliesEditFieldLabel       matlab.ui.control.Label
        NoPliesEditField            matlab.ui.control.NumericEditField
        SpinnerLabel                matlab.ui.control.Label
        Spinner                     matlab.ui.control.Spinner
        PlythicknessEditFieldLabel  matlab.ui.control.Label
        PlythicknessEditField       matlab.ui.control.NumericEditField
        RightPanel                  matlab.ui.container.Panel
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {480, 480};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {446, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 588 480];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {446, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create E1EditFieldLabel
            app.E1EditFieldLabel = uilabel(app.LeftPanel);
            app.E1EditFieldLabel.HorizontalAlignment = 'right';
            app.E1EditFieldLabel.Position = [30 435 25 22];
            app.E1EditFieldLabel.Text = 'E1';

            % Create E1EditField
            app.E1EditField = uieditfield(app.LeftPanel, 'numeric');
            app.E1EditField.Position = [69 435 100 22];

            % Create E2EditFieldLabel
            app.E2EditFieldLabel = uilabel(app.LeftPanel);
            app.E2EditFieldLabel.HorizontalAlignment = 'right';
            app.E2EditFieldLabel.Position = [30 396 25 22];
            app.E2EditFieldLabel.Text = 'E2';

            % Create E2EditField
            app.E2EditField = uieditfield(app.LeftPanel, 'numeric');
            app.E2EditField.Position = [69 396 100 22];

            % Create v12EditFieldLabel
            app.v12EditFieldLabel = uilabel(app.LeftPanel);
            app.v12EditFieldLabel.HorizontalAlignment = 'right';
            app.v12EditFieldLabel.Position = [30 355 25 22];
            app.v12EditFieldLabel.Text = 'v12';

            % Create v12EditField
            app.v12EditField = uieditfield(app.LeftPanel, 'numeric');
            app.v12EditField.Position = [69 355 100 22];

            % Create G12EditFieldLabel
            app.G12EditFieldLabel = uilabel(app.LeftPanel);
            app.G12EditFieldLabel.HorizontalAlignment = 'right';
            app.G12EditFieldLabel.Position = [27 302 28 22];
            app.G12EditFieldLabel.Text = 'G12';

            % Create G12EditField
            app.G12EditField = uieditfield(app.LeftPanel, 'numeric');
            app.G12EditField.Position = [69 302 100 22];

            % Create fxEditFieldLabel
            app.fxEditFieldLabel = uilabel(app.LeftPanel);
            app.fxEditFieldLabel.HorizontalAlignment = 'right';
            app.fxEditFieldLabel.Position = [30 262 25 22];
            app.fxEditFieldLabel.Text = 'fx';

            % Create fxEditField
            app.fxEditField = uieditfield(app.LeftPanel, 'numeric');
            app.fxEditField.Position = [69 262 100 22];

            % Create fyEditFieldLabel
            app.fyEditFieldLabel = uilabel(app.LeftPanel);
            app.fyEditFieldLabel.HorizontalAlignment = 'right';
            app.fyEditFieldLabel.Position = [30 229 25 22];
            app.fyEditFieldLabel.Text = 'fy';

            % Create fyEditField
            app.fyEditField = uieditfield(app.LeftPanel, 'numeric');
            app.fyEditField.Position = [69 229 100 22];

            % Create fxyEditFieldLabel
            app.fxyEditFieldLabel = uilabel(app.LeftPanel);
            app.fxyEditFieldLabel.HorizontalAlignment = 'right';
            app.fxyEditFieldLabel.Position = [30 187 25 22];
            app.fxyEditFieldLabel.Text = 'fxy';

            % Create fxyEditField
            app.fxyEditField = uieditfield(app.LeftPanel, 'numeric');
            app.fxyEditField.Position = [69 187 100 22];

            % Create thetaEditFieldLabel
            app.thetaEditFieldLabel = uilabel(app.LeftPanel);
            app.thetaEditFieldLabel.HorizontalAlignment = 'right';
            app.thetaEditFieldLabel.Position = [23 101 32 22];
            app.thetaEditFieldLabel.Text = 'theta';

            % Create thetaEditField
            app.thetaEditField = uieditfield(app.LeftPanel, 'numeric');
            app.thetaEditField.Position = [69 101 100 22];

            % Create thetaEditField_2Label
            app.thetaEditField_2Label = uilabel(app.LeftPanel);
            app.thetaEditField_2Label.HorizontalAlignment = 'right';
            app.thetaEditField_2Label.Position = [23 68 32 22];
            app.thetaEditField_2Label.Text = 'theta';

            % Create thetaEditField_2
            app.thetaEditField_2 = uieditfield(app.LeftPanel, 'numeric');
            app.thetaEditField_2.Position = [69 68 100 22];

            % Create thetaEditField_3Label
            app.thetaEditField_3Label = uilabel(app.LeftPanel);
            app.thetaEditField_3Label.HorizontalAlignment = 'right';
            app.thetaEditField_3Label.Position = [23 33 32 22];
            app.thetaEditField_3Label.Text = 'theta';

            % Create thetaEditField_3
            app.thetaEditField_3 = uieditfield(app.LeftPanel, 'numeric');
            app.thetaEditField_3.Position = [69 33 100 22];

            % Create NoPliesEditFieldLabel
            app.NoPliesEditFieldLabel = uilabel(app.LeftPanel);
            app.NoPliesEditFieldLabel.HorizontalAlignment = 'right';
            app.NoPliesEditFieldLabel.Position = [3 141 54 22];
            app.NoPliesEditFieldLabel.Text = 'No. Plies';

            % Create NoPliesEditField
            app.NoPliesEditField = uieditfield(app.LeftPanel, 'numeric');
            app.NoPliesEditField.Position = [71 141 100 22];

            % Create SpinnerLabel
            app.SpinnerLabel = uilabel(app.LeftPanel);
            app.SpinnerLabel.HorizontalAlignment = 'right';
            app.SpinnerLabel.Position = [213 33 47 22];
            app.SpinnerLabel.Text = 'Spinner';

            % Create Spinner
            app.Spinner = uispinner(app.LeftPanel);
            app.Spinner.Position = [275 33 100 22];

            % Create PlythicknessEditFieldLabel
            app.PlythicknessEditFieldLabel = uilabel(app.LeftPanel);
            app.PlythicknessEditFieldLabel.HorizontalAlignment = 'right';
            app.PlythicknessEditFieldLabel.Position = [184 435 76 22];
            app.PlythicknessEditFieldLabel.Text = 'Ply thickness';

            % Create PlythicknessEditField
            app.PlythicknessEditField = uieditfield(app.LeftPanel, 'numeric');
            app.PlythicknessEditField.Position = [275 435 100 22];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end