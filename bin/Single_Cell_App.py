#!/bin/python3
"""
Single_Cell_App.py v0.1.0
    Aug. 29, 2018
    Dennis A. Simpson
    This is a functional beta version.  Will work as expected.  Still needs cleaner GUI possibly with larger fonts.
    The Target Search panel output still needs validation.


@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2018
"""
import datetime
import os
import collections
from contextlib import suppress
import dill
import wx
import wx.adv
from wx.lib.wordwrap import wordwrap
import wx.lib.sized_controls as sized_controls
import wx.lib.masked as masked
import wx.lib.intctrl
import subprocess

__author__ = 'Dennis A. Simpson'
__version__ = '0.1.0'
__package__ = 'Volur'
__copyright__ = '(C) 2018'


class IntBoxes:
    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_int_caller(self, name, value=None):
        int_ctrl = wx.lib.intctrl.IntCtrl(self.main_panel, wx.ID_ANY, name=name, value=value, allow_none=1, min=0)
        int_ctrl.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        int_ctrl.SetColors(default_color=wx.BLACK, oob_color=wx.RED)
        return int_ctrl


class DataBoxes:
    """
    Generic ComboBox generator for project.  Ensures all boxes will look the same.
    """

    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_box(self, data_list, name, value=None):
        if value is None:
            value = ""
        my_combo = wx.ComboBox(self.main_panel, wx.ID_DEFAULT, value, choices=data_list, style=wx.CB_SORT, name=name)
        my_combo.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        return my_combo


class MaskedDataBoxes:
    """
    Generic Masked ComboBox generator for project.  Ensures all boxes will look the same.
    """

    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_masked_box(self, data_list, name, value=None):
        my_masked_combo = \
            masked.ComboBox(self.main_panel, wx.ID_DEFAULT, choices=data_list, style=wx.CB_SORT, name=name)
        my_masked_combo.SetCtrlParameters(invalidBackgroundColour="red", defaultValue=value, choiceRequired=True)
        my_masked_combo.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        return my_masked_combo


class ButtonGenerator:
    """
    Generic Button Generator for Project Insuring all Buttons are the same.
    """

    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_buttons(self, label):
        my_button = wx.Button(self.main_panel, wx.ID_ANY, label=label)
        my_button.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        return my_button


class BooleanValidator(wx.Validator):
    def __init__(self, parent, name):
        super(BooleanValidator, self).__init__()
        self.parent = parent
        self.name = name

    def Clone(self):
        return BooleanValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if "/" in text:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} must be a full path statement".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class PathValidator(wx.Validator):
    def __init__(self, parent, name):
        super(PathValidator, self).__init__()
        self.parent = parent
        self.name = name

    def Clone(self):
        """ """
        return PathValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if "/" in text:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} requires full path statements".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class TextValidator(wx.Validator):
    def __init__(self, parent, name):
        super(TextValidator, self).__init__()
        self.parent = parent
        self.name = name

    def Clone(self):
        """ """
        return TextValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if len(text) > 0:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} cannot be Null".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
        return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class NumValidator(wx.lib.intctrl.IntValidator):
    def __init__(self, parent, name):
        wx.lib.intctrl.IntValidator.__init__(self)
        self.parent = parent
        self.name = name

    def Clone(self):
        """ """
        return NumValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if text >= 0:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} takes positive integers only".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()

            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class CommonControls:
    def __init__(self, parent):
        self.parent = parent
        self.default_data = wx.GetTopLevelParent(parent).default_data
        self.dataframe = wx.GetTopLevelParent(parent).dataframe
        self.button_generator = ButtonGenerator(parent)
        self.databoxes = DataBoxes(parent)
        self.my_masked_databoxes = MaskedDataBoxes(parent)
        self.my_int_box = IntBoxes(parent)

    def folder_selector(self, name, tip=None):
        input_ctrl = self.databoxes.my_box(self.dataframe[name], name)
        input_ctrl.SetValidator(PathValidator(self.parent, input_ctrl.GetName()))
        if tip is not None:
            input_ctrl.SetToolTip(tip)
        working_dir_btn = self.button_generator.my_buttons(name)
        working_dir_btn.Bind(wx.EVT_BUTTON, lambda event, temp=working_dir_btn: self.select_folder(event, input_ctrl))

        return working_dir_btn, input_ctrl

    def file_selector(self, name, tip=None):
        input_ctrl = self.databoxes.my_box(self.dataframe[name], name)
        input_ctrl.SetValidator(PathValidator(self.parent, input_ctrl.GetName()))
        if tip is not None:
            input_ctrl.SetToolTip(tip)
        fastq_file_btn = self.button_generator.my_buttons(name)
        fastq_file_btn.Bind(wx.EVT_BUTTON, lambda event, temp=fastq_file_btn: self.select_file(event, input_ctrl))

        return fastq_file_btn, input_ctrl

    def restricted_selector(self, name, value=None, tip=None):
        """
        Define the verbosity level for the logger.
        :return:
        """
        input_ctrl = self.my_masked_databoxes.my_masked_box(self.default_data[name], name, value)
        if tip is not None:
            input_ctrl.SetToolTip(tip)
        label = wx.StaticText(self.parent, wx.ID_ANY, name.strip("--"), style=wx.ALIGN_RIGHT)
        label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        label.SetForegroundColour("blue")
        label.SetMinSize((200, 20))

        return label, input_ctrl

    def text_control(self, name, tip=None):
        try:
            value = self.default_data[name]
        except KeyError:
            value = ""
        input_ctrl = self.databoxes.my_box(self.dataframe[name], name, value)
        input_ctrl.SetValidator(TextValidator(self.parent, input_ctrl.GetName()))
        if tip is not None:
            input_ctrl.SetToolTip(tip)

        label = wx.StaticText(self.parent, wx.ID_DEFAULT, name.strip("--"), style=wx.ALIGN_RIGHT)
        label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        label.SetForegroundColour("blue")
        label.SetMinSize((200, 20))

        return label, input_ctrl

    def int_control(self, name, tip=None):
        input_ctrl = self.my_int_box.my_int_caller(name, self.default_data[name])
        if tip is not None:
            input_ctrl.SetToolTip(tip)

        label = wx.StaticText(self.parent, wx.ID_DEFAULT, name.strip("--"), style=wx.ALIGN_RIGHT)
        label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        # label.SetMinSize((int(self.parent.GetSize()[0] * 0.5), 45))
        label.SetForegroundColour("blue")
        label.SetMinSize((200, 20))

        return label, input_ctrl

    def select_folder(self, event, ctrl):
        dlg = wx.DirDialog(self.parent, "Choose a Folder")
        if dlg.ShowModal() == wx.ID_OK:
            ctrl.SetValue(dlg.GetPath())
            dlg.Destroy()

    def select_file(self, event, ctrl):
        dlg = wx.FileDialog(self.parent, "Choose a file")
        if dlg.ShowModal() == wx.ID_OK:
            ctrl.SetValue(dlg.GetPath())
            dlg.Destroy()

    def int_control_builder(self, name, default_value):
        name = name
        try:
            value = default_value[name]
        except KeyError:
            value = None

        input_ctrl = self.my_int_box.my_int_caller(name, value)

        return input_ctrl


class SingleCellPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(SingleCellPanel, self).__init__(parent, wx.ID_ANY, name="SingleCellPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Single Cell Parameters")

        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list =\
            [my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),
             my_controls.file_selector("--Master_Index_File"),
             my_controls.file_selector("--Index_File"),
             my_controls.restricted_selector("--Verbose", "INFO"),
             my_controls.text_control("--Job_Name"),
             my_controls.int_control("--Spawn"),
             my_controls.restricted_selector("--Species"),
             my_controls.text_control("--Control_Sample"),
             ]

        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
            panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)
            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class SegmentAnalyzerPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(SegmentAnalyzerPanel, self).__init__(parent, wx.ID_ANY, name="SegmentAnalyzerPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Segment Analyzer Parameters")
        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list =\
            [my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),
             my_controls.file_selector("--Segment_File", tip="Full Path to Ginkgo Segment File"),
             my_controls.file_selector("--Target_File", tip="Full Path to Target Bed File"),
             ["Row Skip", ""],
             my_controls.restricted_selector("--Verbose", "INFO"),
             my_controls.text_control("--Job_Name", tip="Used in output.  No spaces or special characters"),
             my_controls.int_control("--Spawn", tip="How many processors/threads to use.  Max should be n-1"),
             my_controls.restricted_selector("--Species"),
             my_controls.text_control("--Cell_Name"),
             my_controls.restricted_selector("--Include_chrY"),
             my_controls.int_control("--Minimum_Family_Size",
                                     tip="Used for Ensemble Gene File. Minimum number of occurrences"),
             ["Row Skip", ""],
             my_controls.restricted_selector("--Breakpoint_Dist_File",
                                             tip="Breakpoint counts by chromosome and type"),
             my_controls.restricted_selector("--Breakpoint_Chrom_Dist_File",
                                             tip="ALL Cells.  Breakpoint Counts by type and chromosome"),
             my_controls.restricted_selector("--Chrom_Ploidy_File", tip="ALL CELLS"),
             my_controls.restricted_selector("--Aberration_Size_File",
                                             tip="ALL CELLS by chromosome and type"),
             my_controls.restricted_selector("--Ensembl_Gene_File",
                                             tip="Reports genes found in targets that intersect"),
             my_controls.restricted_selector("--Target_Intersect_File",
                                             tip="writes data about target intersects"),
             my_controls.restricted_selector("--Breakpoint_Coordinate_File",
                                             tip="Writes files of breakpoint coordinates")
             ]

        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        static_line_sizer = wx.BoxSizer(wx.HORIZONTAL)
        static_line_sizer.Add(wx.StaticLine(self), 0, wx.EXPAND)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            if widget0 == "Row Skip":
                panel_grid_sizer.Add(width=200, height=20, pos=(i, 1), flag=wx.EXPAND)
            else:
                # panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
                panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
                panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)

            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class SamplePermutationPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(SamplePermutationPanel, self).__init__(parent, wx.ID_ANY, name="SamplePermutationPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Sample Permutation Analysis Parameters")
        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list =\
            [my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),
             my_controls.file_selector("--Segment_File", tip="Full Path to Ginkgo Segment File"),
             my_controls.file_selector("--Target_File", tip="Full Path to Target Bed File"),
             ["Row Skip", ""],
             my_controls.restricted_selector("--Verbose", "INFO"),
             my_controls.text_control("--Job_Name", tip="Used in output.  No spaces or special characters"),
             my_controls.int_control("--Spawn", tip="How many processors/threads to use.  Max should be n-1"),
             my_controls.restricted_selector("--Species"),

             my_controls.restricted_selector("--Include_chrY"),
             my_controls.int_control("--Iteration_Count", tip="How iterations to run"),
             my_controls.int_control("--Sample_Group_Size", tip="Size of Sample Permutation Groups."),
             ["Row Skip", ""],
             my_controls.restricted_selector("--Breakpoint_Coordinate_File",
                                             tip="Writes files of breakpoint coordinates."),
             my_controls.int_control("--Prog_Check", tip="How often the program reports its progress")
             ]

        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        static_line_sizer = wx.BoxSizer(wx.HORIZONTAL)
        static_line_sizer.Add(wx.StaticLine(self), 0, wx.EXPAND)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            if widget0 == "Row Skip":
                panel_grid_sizer.Add(width=200, height=20, pos=(i, 1), flag=wx.EXPAND)
            else:
                # panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
                panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
                panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)

            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class SegmentPermutationPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(SegmentPermutationPanel, self).__init__(parent, wx.ID_ANY, name="SegmentPermutationPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Segment Permutation Analysis Parameters")
        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list =\
            [my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),
             my_controls.file_selector("--Segment_File", tip="Full Path to Ginkgo Segment File"),
             my_controls.file_selector("--Target_File", tip="Full Path to Target Bed File"),
             ["Row Skip", ""],
             my_controls.restricted_selector("--Verbose", "INFO"),
             my_controls.text_control("--Job_Name", tip="Used in output.  No spaces or special characters"),
             my_controls.int_control("--Spawn", tip="How many processors/threads to use.  Max should be n-1"),
             my_controls.restricted_selector("--Species"),
             my_controls.text_control("--Cell_Name"),
             my_controls.restricted_selector("--Include_chrY"),
             my_controls.int_control("--Iteration_Count", tip="How iterations to run"),
             my_controls.text_control("--Total_Targeted", tip="Total Observed Target Intersects"),
             my_controls.text_control("--Unique_Targeted", tip="Unique Observed Target Intersects"),
             ["Row Skip", ""],
             my_controls.int_control("--Combine_Segments",
                                     tip="How many segments to combine to simulate breakpoints."),
             my_controls.restricted_selector("--Segment_Permutation_File"),
             my_controls.restricted_selector("--Map_File"),
             ["Row Skip", ""],
             my_controls.int_control("--Prog_Check", tip="How often the program reports its progress")
             ]

        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        static_line_sizer = wx.BoxSizer(wx.HORIZONTAL)
        static_line_sizer.Add(wx.StaticLine(self), 0, wx.EXPAND)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            if widget0 == "Row Skip":
                panel_grid_sizer.Add(width=200, height=20, pos=(i, 1), flag=wx.EXPAND)
            else:
                # panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
                panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
                panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)

            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class WelcomePanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(WelcomePanel, self).__init__(parent, wx.ID_ANY, name="WelcomePanel")
        self.SetBackgroundStyle(wx.BG_STYLE_ERASE)
        self.box = wx.BoxSizer(wx.VERTICAL)

        # m_text = wx.StaticText(self, wx.ID_ANY, "Synthetic Lethal Analysis")
        # m_text.SetFont(wx.Font(20, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.BOLD, faceName='Inconsolata'))
        # m_text.SetForegroundColour("blue")
        # m_text.SetSize(m_text.GetBestSize())
        # self.box.Add(m_text, 0, wx.ALL, 10)
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.add_background)

        self.SetSizerAndFit(self.box)

    def add_background(self, event):
        """
         Add a picture to the background
         """
        dc = wx.ClientDC(self)
        rect = self.GetUpdateRegion().GetBox()
        dc.SetClippingRegion(rect)

        bmp = wx.Bitmap("Asset5.png")
        cliWidth, cliHeight = self.GetClientSize()

        bmpWide = bmp.GetWidth()
        bmpHeight = bmp.GetHeight()
        img = bmp.ConvertToImage()
        bmpWide = int((cliWidth / bmpWide) * bmpWide)
        bmpHeight = int((cliHeight / bmpHeight) * bmpHeight)

        bmp = wx.Bitmap(img.Scale(bmpWide, bmpHeight))

        xPos = (cliWidth - (bmpWide)) / 2
        yPos = (cliHeight - (bmpHeight)) / 2
        dc.DrawBitmap(bmp, xPos, yPos)


class MainFrame(wx.Frame):
    def __init__(self):

        # Setup a frame that is centered on the screen and opens at 75% of full screen.
        window_name = "Single Cell Data Processing GUI  v{}".format(__version__)
        display_x, display_y = wx.DisplaySize()
        scale = 0.75

        wx.Frame.__init__(self, None, wx.ID_ANY, title=window_name, size=(display_x*scale, display_y*scale))
        self.Centre()
        self.build_menu_bar()

        self.dataframe = self.dataframe_build()
        self.default_data = self.default_dict_build()
        print(self.GetSize(), wx.DisplaySize(), wx.version())

        self.welcome_panel = WelcomePanel(self)
        self.segment_analyzer_panel = SegmentAnalyzerPanel(self)
        self.single_cell_panel = SingleCellPanel(self)
        self.segment_permutation_panel = SegmentPermutationPanel(self)
        self.sample_permutation_panel = SamplePermutationPanel(self)

        self.segment_analyzer_panel.Hide()
        self.single_cell_panel.Hide()
        self.segment_permutation_panel.Hide()
        self.sample_permutation_panel.Hide()

        self.panel_dict = \
            {"SegmentAnalyzerPanel": self.segment_analyzer_panel,
             "WelcomePanel": self.welcome_panel,
             "SingleCellPanel": self.single_cell_panel,
             "SegmentPermutationPanel": self.segment_permutation_panel,
             "SamplePermutationPanel": self.sample_permutation_panel
             }

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.panel = self.welcome_panel
        self.main_sizer.Add(self.panel, 1, wx.EXPAND)

        self.SetSizer(self.main_sizer)

    @staticmethod
    def default_dict_build():
        """
        This builds a dictionary with values that are library type and version specific.
        :return:
        """
        version1_default_dict = \
            {"--AnchorStart": 100, "--AnchorStop": 125, "--Minimum_Family_Size": 1, "--Species": ["Human", "Mouse"],
             "--Verbose": ["INFO", "DEBUG", "ERROR"], "--Spawn": 1, "--Include_chrY": ["True", "False"],
             "--Chrom_Ploidy_File": ["True", "False"], "--Breakpoint_Dist_File": ["True", "False"],
             "--Prog_Check": 10000,
             "--Target_Intersect_File": ["True", "False"], "--Aberration_Size_File": ["True", "False"],
             "--Breakpoint_Chrom_Dist_File": ["True", "False"], "--Ensembl_Gene_File": ["True", "False"],
             "--Map_File": ["True", "False"],
             "--Breakpoint_Coordinate_File": ["True", "False"], "--Permutation_Type": ["Segment", "Sample"],
             "--Combine_Segments": 2, "--Segment_Permutation_File": ["True", "False"], "--Iteration_Count": 10000,
             "--Sample_Group_Size": 1
             }
        return version1_default_dict

    @staticmethod
    def dataframe_build():
        dir_path = os.path.dirname(os.path.abspath(__file__))
        pickle_file = "{0}{1}pickles{1}Volur_parameters.pkl".format(dir_path, os.sep)
        try:
            with open(pickle_file, 'rb') as file:
                dataframe_dict = dill.load(file)

        except FileNotFoundError:
            column_names = \
                ["--FASTQ1", "--Index_File", "--Target_File", "--Master_Index_File", "--Working_Folder", "--Verbose",
                 "--Job_Name", "--Spawn", "--Segment_File", "--Include_chrY", "--Cell_Name", "--Combine_Segments",
                 "--Minimum_Family_Size", "--Species", "--Breakpoint_Dist_File", "--Breakpoint_Chrom_Dist_File",
                 "--Chrom_Ploidy_File", "--Aberration_Size_File", "--Ensembl_Gene_File", "--Segment_Permutation_File"
                 "--Target_Intersect_File", "--Write_Map_File", "--Breakpoint_Coordinate_File", "--Iteration_Count",
                 "--Sample_Group_Size", "--Total_Targeted", "--Unique_Targeted"
                 ]

            dataframe_dict = collections.defaultdict(list)
            for key in column_names:
                dataframe_dict[key] = []
        return dataframe_dict

    def build_menu_bar(self):
        menu_bar = wx.MenuBar()

        # Build the "File" menu.
        file_menu = wx.Menu()
        welcome_app = file_menu.Append(wx.ID_ANY, "Home", "Display Welcome Screen")
        save_app = file_menu.Append(wx.ID_SAVE, "&Save\tCtrl-S", "Save Parameter File to Working Folder")
        run_app = file_menu.Append(wx.ID_EXECUTE, "&Run\tCtrl-R", "Save Parameter File and Run")
        file_menu.AppendSeparator()
        exit_app = file_menu.Append(wx.ID_EXIT, "E&xit\tAlt-x", "Close window and exit program.")
        menu_bar.Append(file_menu, "&File")

        # Build the "Tools" Menu
        tools_menu = wx.Menu()
        permutation_submenu = wx.Menu()
        sample_permutation_app = permutation_submenu.Append(wx.ID_ANY, "SamplePermutation",
                                                            "Do random sample permutations")
        segment_permutation_app = permutation_submenu.Append(wx.ID_ANY, "SegmentPermutation",
                                                             "Do random segment permutations")
        tools_menu.Append(wx.ID_ANY, "Permutations", permutation_submenu)
        segment_analyzer_app = tools_menu.Append(wx.ID_ANY, "SegmentAnalyzer", "Analyze Segment Data")
        single_cell_app = tools_menu.Append(wx.ID_ANY, "SingleCell", "")
        menu_bar.Append(tools_menu, "&Tools")

        # Build the "Help" menu
        help_menu = wx.Menu()
        about_item = help_menu.Append(wx.ID_ABOUT)
        menu_bar.Append(help_menu, "&Help")
        self.SetMenuBar(menu_bar)

        # Bind Menu Items to Actions
        self.Bind(wx.EVT_MENU, self.save_parameter_file, save_app)
        self.Bind(wx.EVT_MENU, self.on_run, run_app)
        self.Bind(wx.EVT_MENU, self.on_exit, exit_app)
        self.Bind(wx.EVT_MENU, self.on_about, about_item)

        # Panel switching
        self.Bind(wx.EVT_MENU, lambda event, temp=segment_analyzer_app:
                  self.switch_panels(event, self.segment_analyzer_panel.GetName()), segment_analyzer_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=welcome_app:
                  self.switch_panels(event, self.welcome_panel.GetName()), welcome_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=sample_permutation_app:
                  self.switch_panels(event, self.sample_permutation_panel.GetName()), sample_permutation_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=segment_permutation_app:
                  self.switch_panels(event, self.segment_permutation_panel.GetName()), segment_permutation_app)

    def switch_panels(self, event, switch_id):
        self.main_sizer.Detach(self.panel)
        self.panel.Hide()
        self.panel = self.panel_dict[switch_id]
        self.main_sizer.Add(self.panel, 1, wx.EXPAND)
        self.panel.Show()
        self.panel.Fit()
        self.Layout()

    def on_exit(self, event):
        """Close the frame, terminating the application."""
        self.Close(True)

    def on_about(self, event):
        about_info = wx.adv.AboutDialogInfo()

        about_info.SetName("Single Cell App")
        about_info.SetCopyright(__copyright__)

        about_info.SetDescription(wordwrap(
            "GUI for Single Cell Pipeline.\nVersion: {}".format(__version__), 350, wx.ClientDC(self)))
        about_info.SetDevelopers([__author__])
        # about_info.License = wordwrap("Completely and totally open source!", 500, wx.ClientDC(self))
        wx.adv.AboutBox(about_info)

    def on_run(self, event):
        shellfile_name = self.save_parameter_file("Run")
        print(shellfile_name)
        if shellfile_name is not None:
            subprocess.run([shellfile_name], shell=True)

    def save_parameter_file(self, event):
        def outstring_build():
            file_body = "{}{}{}\n".format(do_permutations, segment_analyzer, single_cell)
            working_folder = ""
            job = ""

            for i in range(len(panel.control_dict)):
                d = panel.control_dict[i][1]

                # Validate data.  If validation fails this will take us back to the form.
                with suppress(AttributeError):
                    if not d.GetValidator().Validate(d):
                        return job, working_folder, file_body, False

                try:
                    ctrl_name = d.GetName()
                    ctrl_value = str(d.GetValue()).strip()
                except AttributeError:
                    ctrl_value = ctrl_name = None

                # Linux does not add the trailing slash to folders.  I find that confusing.

                if ctrl_name == "--Working_Folder":
                    working_folder = ctrl_value
                    ctrl_value += "/"
                if ctrl_name == "--Job_Name":
                    job = ctrl_value
                if ctrl_name is None:
                    file_body += "\n"
                else:
                    file_body += "{}\t{}\n".format(ctrl_name, ctrl_value)
                    self.dataframe[ctrl_name].append(d.GetValue())
                    self.dataframe[ctrl_name] = list(set(self.dataframe[ctrl_name]))

            return job, working_folder, file_body, True
        segment_analyzer = "--Segment_Analyzer\tFalse\n"
        single_cell = "--Single_Cell\tFalse\n"
        do_permutations = "--Permutation_Analysis\tFalse\n"
        panel_parameters = ""

        if self.segment_analyzer_panel.IsShown():
            submodule = "Segment_Analyzer"
            panel = self.segment_analyzer_panel
            segment_analyzer = "--Segment_Analyzer\tTrue\n"

        elif self.single_cell_panel.IsShown():
            submodule = "Single_Cell"
            panel = self.single_cell_panel
            single_cell = "--Single_Cell\tTrue\n"

        elif self.sample_permutation_panel.IsShown():
            submodule = "Sample Permutation_Analysis"
            panel = self.sample_permutation_panel
            do_permutations = "--Permutation_Analysis\tTrue\n" \
                              "--Permutation_Type\tSample\n"
            panel_parameters = "--Combine_Segments\t1\n" \
                               "--Map_File\tFalse\n" \
                               "--Segment_Permutation_File\tFalse\n" \
                               "--Target_Intersect_File\tFalse\n"

        elif self.segment_permutation_panel.IsShown():
            submodule = "Segment Permutation_Analysis"
            panel = self.segment_permutation_panel
            do_permutations = "--Permutation_Analysis\tTrue\n" \
                              "--Permutation_Type\tSegment\n"

        job_name, working_dir, parameter_body, validation_pass = outstring_build()

        # If the validation fails take the user back to the panel so they can correct the error(s)
        if not validation_pass:
            return

        dir_path = os.path.dirname(os.path.abspath(__file__))
        date = datetime.datetime.now().strftime("%Y%m%d")
        outfile_name = "{}{}run_{}_{}.sh".format(working_dir, os.sep, job_name, date)
        shebang = "#!/bin/bash\n" \
                  "#Parameter file to run Volur Single Cell module {}\n" \
                  "#File generated {}\n\n".format(submodule, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        cmd = "python3 {0}{1}Volur.py --options_file {2}\nexit\n\n".format(dir_path, os.sep, outfile_name)

        # Save the current choices in the pickles folder.
        pickle_file = "{0}{1}pickles{1}Volur_parameters.pkl".format(dir_path, os.sep)
        with open(pickle_file, 'wb') as file:
            dill.dump(self.dataframe, file, protocol=-1)

        # Write the parameter file.
        outstring = "{}{}{}{}".format(shebang, cmd, parameter_body, panel_parameters)
        parameter_file = \
            open(outfile_name, 'w')
        parameter_file.write(outstring)
        parameter_file.close()

        if event == "Run" and validation_pass:
            return outfile_name


def main():
    app = wx.App()
    frame = MainFrame()
    frame.Show()
    app.MainLoop()
    # wx.lib.inspection.InspectionTool().Show()


if __name__ == '__main__':
    main()
