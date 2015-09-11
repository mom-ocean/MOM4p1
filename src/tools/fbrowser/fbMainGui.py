#!/usr/bin/env python

import sys, os, math, Tkinter, tkFileDialog, types, tkFont, Pmw, Tix
from fbDatabaseFrame import fbDatabaseFrame
from fbExpBrowser import fbExpBrowser
from fbFigBrowser import fbFigBrowser

__FILE__ = sys._getframe().f_code.co_filename

def makeResizable(widget):
    for i in range(0, widget.grid_size()[0]):
        widget.columnconfigure(i, weight=1)
    for i in range(0, widget.grid_size()[1]):
        widget.rowconfigure(i, weight=1)

class fbMainGui:

    def __init__(self, root=None, title=None, initExpDirs=None):
        """
        Constructor
        root: the parent container widget (e.g. Toplevel, Frame, etc.)
        title: the title of this window
        initExpDirs: list of directories to be loaded into the Selected Experiments list
        """
        if root:
            self.root = root
        else:
            self.root = Tix.Tk()

        self.root.geometry("1000x600+0+0")

        if title:
            self.root.title(title)
        else:
            self.root.title("Figure Browser")

        if initExpDirs:
            self.initExpDirs = initExpDirs
        else:
            self.initExpDirs = []

        #get the available space on the screen
        testWin = Tix.Tk()
        testWin.geometry("100000x100000+0+0")
        testWin.update()
        #self.scr_width = testWin.winfo_width()
        self.scr_width = testWin.winfo_screenwidth()
        self.scr_height = testWin.winfo_height()
        testWin.destroy()


        self.fontFamily = "helvetica"
        self.baseSize = 12

        self.baseSizeVar = Tkinter.IntVar()
        self.baseSizeVar.set(12)

        user = os.getenv("USER")
        self.initBrowseDir = "/home/" + user + "/"
        self.configPath = "/home/" + user + "/.fbconfig"

        self.archiveWarningVar = Tkinter.BooleanVar()
        self.archiveWarningVar.set(True)

        self.gvApp = "gv"
        self.gridFigures = True
        self.forceOrder = False
        self.gvOrient = "-orientation=portrait"
        self.gvAntiAlias = "-noantialias"
        #self.orient = "-portrait"
        #self.antiAlias = "-noantialias"
        self.gvScale = 0
        self.mgvOrient = "-portrait"
        self.mgvAntiAlias = "-noantialias"
        self.mgvAutoHeight = "-noautoheight"
        self.mgvAutoWidth = "-noautowidth"
        self.mgvMagStep = 0

        if os.path.isfile(self.configPath):
            file = open(self.configPath, "r")
            lines = file.readlines()
            options = {}
            for line in lines:
                words = line.split("=")

                key = words[0].lower().strip()
                val = words[1].strip().strip("\'\"\n")
                options[key] = val

            if "app" in options.keys():
                val = options["app"]
                app = val.lower()
                if app == "gv" or app == "mgv":
                    self.gvApp = app
                else:
                    print "'app' must be either 'gv' or 'mgv'."

            for key, val in options.iteritems():

                if key == "fontfamily":
                    if val in tkFont.families():
                        self.fontFamily = val

                elif key == "fontsize":
                    self.baseSize = int(val)
                    self.baseSizeVar.set(int(val))

                elif key == "browsedir":
                    browseDir = val

                    if browseDir.lower() == "pwd":
                        self.initBrowseDir = os.getcwd()
                    elif os.path.exists(browseDir):
                        self.initBrowseDir = browseDir

                    if not self.initBrowseDir.endswith("/"):
                        self.initBrowseDir = self.initBrowseDir + "/"

                elif key == "archivewarning":
                    val = val.lower()
                    if val in ["t", "f"]:
                        if val == "t":
                            self.archiveWarningVar.set(True)
                        else:
                            self.archiveWarningVar.set(False)

                elif key == "app":
                    pass

                elif key == "gridfigures":
                    val = val.lower()
                    if val in ["t", "f"]:
                        if val == "t":
                            val = True
                        else:
                            val = False
                        self.gridFigures = val

                elif key == "forceorder":
                    val = val.lower()
                    if val in ["t", "f"]:
                        if val == "t":
                            val = True
                        else:
                            val = False
                        self.forceOrder = val

                elif key == "orient":
                    val = val.lower()
                    if val in ["portrait", "landscape", "upsidedown", "seascape"]:
                        self.gvOrient = "-orientation=" + val
                        self.mgvOrient = "-" + val

                elif key == "antialias":
                    val = val.lower()
                    if val in ["t", "f"]:
                        if val == "t":
                            val = "antialias"
                        else:
                            val = "noantialias"
                        self.gvAntiAlias = "-" + val
                        self.mgvAntiAlias = "-" + val

                elif key == "scale":
                    self.gvScale = int(val)

                elif key == "magstep":
                    self.mgvMagStep = int(val)

                elif key == "autoheight":
                    val = val.lower()
                    if val in ["t", "f"]:
                        if val == "t":
                            self.mgvAutoHeight = "-autoheight"
                        else:
                            self.mgvAutoHeight = "-noautoheight"

                elif key == "autowidth":
                    val = val.lower()
                    if val in ["t", "f"]:
                        if val == "t":
                            self.mgvAutoWidth = "-autowidth"
                        else:
                            self.mgvAutoWidth = "-noautowidth"

                else:
                    print "'" + key + "' is not a recognized option name."

            file.close()
        else:
            file = open(self.configPath, "w")
            file.write("fontfamily = " + self.fontFamily + "\n")
            file.write("fontsize = " + str(self.baseSize) + "\n")
            file.write("browsedir = " + self.initBrowseDir + "\n")
            if self.archiveWarningVar.get():
                file.write("archivewarning = t\n")
            else:
                file.write("archivewarning = f\n")
            file.write("app = " + self.gvApp + "\n")
            file.write("orient = " + self.gvOrient[1:] + "\n")
            if self.gvAntiAlias == "-antialias":
                file.write("antialias = t\n")
            else: #elif self.gvAntiAlias == "-noantialias":
                file.write("antialias = f\n")
            file.write("scale = " + str(self.gvScale) + "\n")
            file.close()

        self.fontSizes = [self.baseSize, self.baseSize+6, self.baseSize-2, self.baseSize, self.baseSize+4]

        self.baseFont       = tkFont.Font(family=self.fontFamily, size=self.fontSizes[0], weight="bold")
        self.messageBarFont = tkFont.Font(family=self.fontFamily, size=self.fontSizes[1], weight="normal")
        self.entryFieldFont = tkFont.Font(family=self.fontFamily, size=self.fontSizes[2], weight="normal")
        self.optLabelFont   = tkFont.Font(family=self.fontFamily, size=self.fontSizes[3], weight="bold",
                                                                                              underline=1)
        self.gvButtonFont   = tkFont.Font(family=self.fontFamily, size=self.fontSizes[4], weight="bold")

        self.root.option_add("*font", self.baseFont)

        self.fonts = [self.baseFont, self.messageBarFont, self.entryFieldFont,
                      self.optLabelFont, self.gvButtonFont]

        self.defaultFonts = [ {"family":self.fontFamily, "size":self.fontSizes[0],
                               "weight":"bold",   "underline":0},
                              {"family":self.fontFamily, "size":self.fontSizes[1],
                               "weight":"normal", "underline":0},
                              {"family":self.fontFamily, "size":self.fontSizes[2],
                               "weight":"normal", "underline":0},
                              {"family":self.fontFamily, "size":self.fontSizes[3],
                               "weight":"bold",   "underline":1},
                              {"family":self.fontFamily, "size":self.fontSizes[4],
                               "weight":"bold",   "underline":0} ]


        self.SetupInstructions()
        self.SetupAboutDialog()
        self.SetupChangeLog()

        self.menuBar = Pmw.MainMenuBar(self.root)
        self.root.configure(menu = self.menuBar)

        self.menuBar.addmenu("Options", "")
        self.menuBar.addmenuitem("Options", "command", label="Font", command=self.ShowFontOptions)
        self.menuBar.addmenuitem("Options", "command", label="GhostView", command=self.ShowGVOptions)
        self.menuBar.addmenuitem("Options", "command", label="Misc", command=self.ShowMiscOptions)
        self.menuBar.addmenuitem("Options", "separator")
        self.menuBar.addmenuitem("Options", "command", label="Save Options to .fbconfig",
                                             command=self.SaveOptions)

        self.menuBar.addmenu("Change Log", "")
        self.menuBar.addmenuitem("Change Log", "command", label="Change Log", command=self.ShowChangeLog)

        self.menuBar.addmenu("Help", "")
        self.menuBar.addmenuitem("Help", "command", label="Instructions", command=self.ShowInstructions)
        self.menuBar.addmenuitem("Help", "command", label="About", command=self.ShowAboutDialog)

	self.menuBar.add_command(label="Quit", command=root.quit)

        self.fontOptionsVisible = False
        self.SetupFontOptions()
        self.fontOptionsWindow.protocol("WM_DELETE_WINDOW", self.HideFontOptions)


        self.gvOptionsVisible = False
        self.SetupGVOptions()
        self.gvOptionsWindow.protocol("WM_DELETE_WINDOW", self.HideGVOptions)

        self.miscOptionsVisible = False
        self.SetupMiscOptions()
        self.miscOptionsWindow.protocol("WM_DELETE_WINDOW", self.HideMiscOptions)

        self.SetupArchiveWarning()


        self.LR_panes = Pmw.PanedWidget(self.root,
                                        #hull_borderwidth=1,
                                        hull_relief='sunken',
                                        orient='horizontal',
                                        handlesize=0,
                                        separatorrelief='raised',
                                        separatorthickness=8)
        self.leftPane = self.LR_panes.add("leftPane", size=0.25)
        self.rightPane = self.LR_panes.add("rightPane", size=0.75)

        self.TB_panes = Pmw.PanedWidget(self.rightPane,
                                        hull_borderwidth=1,
                                        hull_relief='sunken',
                                        orient='vertical',
                                        handlesize=0,
                                        separatorrelief='raised',
                                        separatorthickness=8)
        self.topRightPane = self.TB_panes.add("topRightPane", size=0.3)
        self.bottomRightPane = self.TB_panes.add("bottomRightPane", size=0.7)

        self.dbFrameSelected = Tkinter.BooleanVar()
        self.dbFrameSelected.set(False)
        self.dbFrameSelected.trace_variable('w', self.ShowFigDir)

        self.dbFrameConfirmed = Tkinter.BooleanVar()
        self.dbFrameConfirmed.set(False)
        self.dbFrameConfirmed.trace_variable('w', self.AddFigDir)

        self.expBrowserModified = Tkinter.BooleanVar()
        self.expBrowserModified.set(False)
        self.expBrowserModified.trace_variable('w', self.ShowFileTree)

        self.figBrowserModified = Tkinter.BooleanVar()
        self.figBrowserModified.set(False)
        self.figBrowserModified.trace_variable('w', self.LaunchGhostView)

        self.messageBar = Pmw.MessageBar(self.root,
                                         silent=1,
                                         entry_font = self.messageBarFont)


        self.dbFrame = fbDatabaseFrame(self.leftPane, self.dbFrameSelected, self.dbFrameConfirmed,
                                       self.messageBar)
        self.dbFrame.grid(sticky="nsew")

        self.expBrowser = fbExpBrowser(self.topRightPane, self.expBrowserModified, self.messageBar,
                                       self.entryFieldFont, self.initBrowseDir, self.initExpDirs)
        self.expBrowser.grid(sticky="nsew")

        self.figBrowser = fbFigBrowser(self.bottomRightPane, self.figBrowserModified, self.messageBar,
                                       self.gvButtonFont, self.entryFieldFont)
        self.figBrowser.grid(sticky="nsew")

        self.LR_panes.grid(row=0, rowspan=10, sticky="nsew")
        self.leftPane.grid(sticky="nsew")
        self.rightPane.grid(sticky="nsew")

        self.TB_panes.grid(sticky="nsew")
        self.topRightPane.grid(sticky="nsew")
        self.bottomRightPane.grid(sticky="nsew")

        self.messageBar.grid(row=10, rowspan=1, sticky='nsew')


        makeResizable(self.root)
        makeResizable(self.leftPane)
        makeResizable(self.rightPane)
        makeResizable(self.topRightPane)
        makeResizable(self.bottomRightPane)
        makeResizable(self.fontOptionsWindow)
        makeResizable(self.gvOptionsWindow)

        self.root.rowconfigure(self.root.grid_size()[1]-1, weight=0)

        self.root.update()


    def SetupFontOptions(self):
        self.fontOptionsWindow = Tkinter.Toplevel()
        self.fontOptionsWindow.withdraw()
        self.fontOptionsWindow.title("Font Settings")
        self.fontOptionsWindow.transient(self.root)

        baseFont_const     = tkFont.Font(family=self.fontFamily, size=self.fontSizes[0], weight="bold")
        optLabelFont_const = tkFont.Font(family=self.fontFamily, size=self.fontSizes[3], weight="bold",
                                                                                            underline=1)

        self.FontSizeLabel = Tkinter.Label(self.fontOptionsWindow,
                                           text = "Font Size",
                                           font = optLabelFont_const)
        self.FontUpButton = Tkinter.Button(self.fontOptionsWindow,
                                           text = "Up",
                                           font = baseFont_const,
                                           command = lambda: self.ChangeFontSizes(1))
        self.FontDownButton = Tkinter.Button(self.fontOptionsWindow,
                                             text = "Down",
                                             font = baseFont_const,
                                             command = lambda: self.ChangeFontSizes(-1))

        self.FontValueLabel = Tkinter.Label(self.fontOptionsWindow,
                                            textvariable=self.baseSizeVar,
                                            bg="#d8d8d8",
                                            font = self.baseFont,
                                            relief='ridge')

        self.FontFamilyLabel = Tkinter.Label(self.fontOptionsWindow,
                                             text = "Font Family",
                                             font = optLabelFont_const)

        self.FontFamilyListBox = Pmw.ScrolledListBox(self.fontOptionsWindow,
                                                     #listbox_height=5,
                                                     dblclickcommand=self.ChangeFontFamily,
                                                     listbox_font = baseFont_const,
                                                     listbox_exportselection=0)

        families = list(tkFont.families())
        families.sort()
        cur = 0
        for i in range(len(families)):
            f = families[i]
            self.FontFamilyListBox.insert('end', f)
            if f == self.fontFamily:
                cur = i
        self.FontFamilyListBox.see(cur)
        self.FontFamilyListBox.selection_set(cur)


        self.ResetButton = Tkinter.Button(self.fontOptionsWindow,
                                          text = "Reset Fonts",
                                          font = baseFont_const,
                                          command = self.ResetFonts)


        self.FontSizeLabel.grid(row=0, column=0, sticky='ew')
        self.FontUpButton.grid(row=0, column=1, sticky='ew')
        self.FontValueLabel.grid(row=0, column=2, sticky='nsew')
        self.FontDownButton.grid(row=0, column=3, sticky='ew')

        self.FontFamilyLabel.grid(row=1, column=0, sticky='ew')
        self.FontFamilyListBox.grid(row=1, rowspan=5, column=1, columnspan=3, sticky='nsew')

        self.ResetButton.grid(row=5, column=0, sticky='sew')

        self.root.winfo_toplevel().update_idletasks()
        self.root.update()
        width = self.fontOptionsWindow.winfo_reqwidth()
        height = self.fontOptionsWindow.winfo_reqheight()
        self.fontOptionsWindow.geometry(str(width) + "x" + str(height) + "+"
                                  + str(0) + "+"
                                  + str(self.menuBar.winfo_reqheight()))

    def ChangeFontFamily(self, *args):
        family = self.FontFamilyListBox.getvalue()[0]
        self.fontFamily = family
        for f in self.fonts:
            f.configure(family=self.fontFamily)
        self.root.update()

    def ChangeFontSizes(self, n):
        if n < 0 and self.baseFont.cget("size") <= 4:
            return
        elif n > 0 and self.baseFont.cget("size") >= 128:
            return

        self.baseSizeVar.set(self.fontSizes[0] + n)
        for i in range(len(self.fonts)):
            self.fontSizes[i] += n
            self.fonts[i].configure(size = self.fontSizes[i])
        self.root.update()

    def ResetFonts(self, *args):
        self.fontFamily = self.defaultFonts[0]["family"]
        self.baseSizeVar.set(self.defaultFonts[0]["size"])
        for i in range(len(self.fonts)):
            defaults = self.defaultFonts[i]
            self.fonts[i].configure(family = defaults["family"],
                                    size = defaults["size"],
                                    weight = defaults["weight"],
                                    underline = defaults["underline"])
            self.fontSizes[i] = defaults["size"]
        self.root.update()


    def HideFontOptions(self, *args):
        self.fontOptionsWindow.withdraw()
        self.fontOptionsVisible = False

    def ShowFontOptions(self, *args):
        if not self.fontOptionsVisible:
            self.fontOptionsWindow.deiconify()
        self.root.winfo_toplevel().update_idletasks()


    def ChangeApp(self, *args):
        self.opt_frame.grid_remove()

        curApp = self.appVar.get()
        self.opt_frame = self.both_frames[curApp]
        self.opt_options = self.both_options[curApp]
        self.opt_controlVars = self.both_controlVars[curApp]
        self.opt_prefixes = self.both_prefixes[curApp]

        self.opt_frame.grid()

        self.root.winfo_toplevel().update_idletasks()
        self.root.update()
        width = self.gvOptionsWindow.winfo_reqwidth()
        height = self.gvOptionsWindow.winfo_reqheight()

        x = self.gvOptionsWindow.winfo_x()
        y = self.gvOptionsWindow.winfo_y()
        self.gvOptionsWindow.geometry(str(width) + "x" + str(height) + "+"
                                   + str(x) + "+" + str(y))

    def SetupGVOptions(self):

        self.gvOptionsWindow = Tkinter.Toplevel()
        self.gvOptionsWindow.withdraw()
        self.gvOptionsWindow.title("GhostView Settings")
        self.gvOptionsWindow.transient(self.root)

        #choice between "gv" and "mgv"
        self.appVar = Tkinter.StringVar()
        self.appVar.set(self.gvApp)
        self.appVar.trace_variable("w", self.ChangeApp)

        self.apps = ["gv", "mgv"]

        radio = Tkinter.Radiobutton(self.gvOptionsWindow,
                                    text="gv",
                                    variable=self.appVar,
                                    value="gv")
        radio.grid(row=0, column=4, columnspan=1, sticky='e')
        radio = Tkinter.Radiobutton(self.gvOptionsWindow,
                                    text="mgv",
                                    variable=self.appVar,
                                    value="mgv")
        radio.grid(row=0, column=5, columnspan=1, sticky='w')

        self.gridFiguresVar = Tkinter.BooleanVar()
        self.gridFiguresVar.set(self.gridFigures)

        check = Tkinter.Checkbutton(self.gvOptionsWindow,
                                    text="Display figures in grid layout on screen",
                                    variable=self.gridFiguresVar)
        check.grid(row=1, column=2, columnspan=8, sticky='w')

        self.forceOrderVar = Tkinter.BooleanVar()
        self.forceOrderVar.set(self.forceOrder)

        check = Tkinter.Checkbutton(self.gvOptionsWindow,
                                    text="Force figures to appear in order of selection",
                                    variable=self.forceOrderVar)
        check.grid(row=2, column=2, columnspan=8, sticky='w')


        self.both_frames = {}
        self.both_frames["gv"] = Tkinter.Frame(self.gvOptionsWindow)
        self.both_frames["mgv"] = Tkinter.Frame(self.gvOptionsWindow)

        self.both_options = {}
        self.both_options["gv"] = [

          #Name             Widget Type     Var Type        Possible values               Default value
          #                                            (Full-text    , command-line)
          ("Orientation",   "radiobutton",  "string", [("Portrait"   , "-orientation=portrait"),
                                                       ("Landscape"  , "-orientation=landscape"),
                                                       ("Upside-Down", "-orientation=upside-down"),
                                                       ("Seascape"   , "-orientation=seascape")],     self.gvOrient),

          ("Antialias",     "radiobutton",  "string", [("On"         , "-antialias"),
                                                       ("Off"        , "-noantialias")],  self.gvAntiAlias),

          ("Scale",         "scale",        "int",    (-5, 5),        "-scale=",           self.gvScale)

        ]

        self.both_options["mgv"] = [

          #Name           Widget Type     Var Type        Possible values               Default value
          #                                            (Full-text    , command-line)
          ("Orientation", "radiobutton",  "string", [("Portrait"   , "-portrait"),
                                                     ("Landscape"  , "-landscape"),
                                                     ("Upside-Down", "-upsidedown"),
                                                     ("Seascape"   , "-seascape")],     self.mgvOrient),

          ("Antialias",   "radiobutton",  "string", [("On"         , "-antialias"),
                                                     ("Off"        , "-noantialias")],  self.mgvAntiAlias),

          ("AutoHeight",  "radiobutton",  "string", [("On"         , "-autoheight"),
                                                     ("Off"        , "-noautoheight")], self.mgvAutoHeight),

          ("AutoWidth",   "radiobutton",  "string", [("On"         , "-autowidth"),
                                                     ("Off"        , "-noautowidth")],  self.mgvAutoWidth),

          ("Mag. Step",   "scale",        "int",     (-5, 5),        "-magstep",        self.mgvMagStep)

        ]

        self.both_controlVars = {}
        self.both_controlVars["gv"] = {}
        self.both_controlVars["mgv"] = {}

        self.both_prefixes = {}
        self.both_prefixes["gv"] = {}
        self.both_prefixes["mgv"] = {}


        for appName in self.apps:
            row=1
            col=0
            options = self.both_options[appName]
            controlVars = self.both_controlVars[appName]
            prefixes = self.both_prefixes[appName]
            frame = self.both_frames[appName]
            for entry in options:
                col=0

                myName = entry[0]
                myWidgetType = entry[1]
                myVarType = entry[2]
                myDefault = entry[-1]

                if myVarType == "int":
                    controlVars[myName] = Tkinter.IntVar()
                elif myVarType == "string":
                    controlVars[myName] = Tkinter.StringVar()
                controlVars[myName].set(myDefault)

                label = Tkinter.Label(frame,
                                      text=myName,
                                      font = self.optLabelFont)

                if myWidgetType == "scale":
                    label.grid(row=row, column=0, sticky='s')
                else:
                    label.grid(row=row, column=0)
                col = col+1

                if myWidgetType == "radiobutton":
                    myPossible = entry[3]
                    for item in myPossible:
                        radio = Tkinter.Radiobutton(frame,
                                                    text=item[0],
                                                    variable=controlVars[myName],
                                                    value=item[1])
                        radio.grid(row=row, column=col, sticky='w')
                        col = col+1
                        prefixes[myName] = ""

                elif myWidgetType == "scale":
                    myRange = entry[3]
                    scale = Tkinter.Scale(frame,
                                          variable=controlVars[myName],
                                          orient='horizontal',
                                          from_=myRange[0],
                                          to=myRange[1])
                    scale.grid(row=row, column=col, columnspan=5, sticky='ew')
		    # HACK ALERT - cmr
		    # Ghostview's command line has apparently changed since this was written
		    # instead of -scale x it wants -scale=x where X is the scaling value. 
		    # mgv, fortunately, has not changed (fortune is a relative term, though.
		    # since the two behave differently it needs this hack to get the arguments correct.
		    if entry[4] == "-scale=":
                        prefixes[myName] = entry[4]
		    else:
			prefixes[myName] = entry[4] + " "

                row = row+1

            makeResizable(frame)

        curApp = self.appVar.get()
        self.opt_frame = self.both_frames[curApp]
        self.opt_options = self.both_options[curApp]
        self.opt_controlVars = self.both_controlVars[curApp]
        self.opt_prefixes = self.both_prefixes[curApp]

        for f in self.both_frames.keys():
            self.both_frames[f].grid(row=3, column=0, columnspan=10)
            self.both_frames[f].grid_remove()

        self.opt_frame.grid()
        makeResizable(self.gvOptionsWindow)

        self.root.winfo_toplevel().update_idletasks()
        self.root.update()
        width = self.gvOptionsWindow.winfo_reqwidth()
        height = self.gvOptionsWindow.winfo_reqheight()
        self.gvOptionsWindow.geometry(str(width) + "x" + str(height) + "+"
                                  + str(self.fontOptionsWindow.winfo_reqwidth() + 10) + "+"
                                  + str(self.menuBar.winfo_reqheight()))

    def HideGVOptions(self, *args):
        self.gvOptionsWindow.withdraw()
        self.gvOptionsVisible = False

    def ShowGVOptions(self, *args):
        if not self.gvOptionsVisible:
            self.gvOptionsWindow.deiconify()
        self.root.winfo_toplevel().update_idletasks()

    def SetupMiscOptions(self, *args):
        self.miscOptionsWindow = Tkinter.Toplevel()
        self.miscOptionsWindow.withdraw()
        self.miscOptionsWindow.title("Misc Settings")
        self.miscOptionsWindow.transient(self.root)

        self.warning_Checkbutton = Tkinter.Checkbutton(self.miscOptionsWindow,
                                                       text="Warn when listing from /archive",
                                                       variable=self.archiveWarningVar)
        self.warning_Checkbutton.grid(row=0, rowspan=2, column=0, columnspan=2, sticky='ew')

        makeResizable(self.miscOptionsWindow)

    def HideMiscOptions(self, *args):
        self.miscOptionsWindow.withdraw()
        self.miscOptionsVisible = False

    def ShowMiscOptions(self, *args):
        if not self.miscOptionsVisible:
            self.miscOptionsWindow.deiconify()
        self.root.winfo_toplevel().update_idletasks()

    def SaveOptions(self, *args):
        file = open(self.configPath, "w")

        file.write("fontfamily = " + self.fontFamily + "\n")

        file.write("fontsize = " + str(self.baseSizeVar.get()) + "\n")

        file.write("browsedir = " + self.initBrowseDir + "\n")

        if self.archiveWarningVar.get():
            file.write("archivewarning = t\n")
        else:
            file.write("archivewarning = f\n")


        app = self.appVar.get()
        file.write("app = " + app + "\n")

        if self.gridFiguresVar.get():
            file.write("gridfigures = t\n")
        else:
            file.write("gridfigures = f\n")

        if self.forceOrderVar.get():
            file.write("forceorder = t\n")
        else:
            file.write("forceorder = f\n")


        controlVars = self.both_controlVars[app]

        file.write("orient = " + controlVars["Orientation"].get()[1:] + "\n")

        if controlVars["Antialias"].get() == "-antialias":
            file.write("antialias = t\n")
        else: #elif controlVars["Antialias"].get() == "-noantialias":
            file.write("antialias = f\n")

        file.write("scale = " + str((self.both_controlVars["gv"])["Scale"].get()) + "\n")
        file.write("magstep = " + str((self.both_controlVars["mgv"])["Mag. Step"].get()) + "\n")

        controlVars = self.both_controlVars["mgv"]

        if controlVars["AutoWidth"].get() == "-autowidth":
            file.write("autowidth = t\n")
        else: #elif controlVars["AutoWidth"].get() == "-noautowidth":
            file.write("autowidth = f\n")

        if controlVars["AutoHeight"].get() == "-autoheight":
            file.write("autoheight = t\n")
        else: #elif controlVars["AutoHeight"].get() == "-noautoheight":
            file.write("autoheight = f\n")

        file.close()

    def SetupInstructions(self, *args):
        self.instructWindow = Pmw.TextDialog(self.root,
                                             title = "Instructions",
                                             text_wrap='word',
                                             #text_width=85,
                                             )
        self.instructWindow.withdraw()
        self.instructWindow.transient(self.root)

        self.instructWindow.tag_config("level1_text", lmargin1=0, lmargin2=0)
        self.instructWindow.tag_config("level1_listhead", lmargin1="0.5c", lmargin2="0.5c")
        self.instructWindow.tag_config("level2_text", lmargin1="1c", lmargin2="1c")
        self.instructWindow.tag_config("level2_listhead", lmargin1="1.5c", lmargin2="1.5c")
        self.instructWindow.tag_config("level3_text", lmargin1="2c", lmargin2="2c")

        self.instructWindow.insert('end',
            ("The Figure Browser is an application for easily browsing, selecting, and viewing" +
            " figures contained in PostScript and PDF files.\n\n" +
            "The Figure Browser's window is divided into 3 panes:\n\n"),
        "level1_text")

        self.instructWindow.insert('end', "1) ", "level1_listhead")
        self.instructWindow.insert('end',
            ("The LEFT pane is called the DATABASE BROWSER. It lists all the experiments that" +
            " you can view.\n"),
        "level2_text")
        self.instructWindow.insert('end', "2) ", "level1_listhead")
        self.instructWindow.insert('end',
            ("The UPPER-RIGHT pane is called the EXPERIMENT SELECTOR. Here you can view a" +
            " list of the experiments you've already selected, and add or delete experiments" +
            " from the list.\n"),
        "level2_text")
        self.instructWindow.insert('end', "3) ", "level1_listhead")
        self.instructWindow.insert('end',
            ("The LOWER-RIGHT pane is called the FIGURE SELECTOR. Here you can select a figure" +
            " to be opened in GhostView.\n\n"),
        "level2_text")

        self.instructWindow.insert('end',
            "Using the application involves 3 basic steps:\n\n",
        "level1_text")

        self.instructWindow.insert('end', "1) ", "level1_listhead")
        self.instructWindow.insert('end',
            ("Select the experiments that you wish to view; the EXPERIMENT SELECTOR lists the" +
            " experiments you have already selected. You can add a new experiment to the list" +
            " in one of 4 ways:\n\n"),
        "level2_text")

        self.instructWindow.insert('end', "a) ", "level2_listhead")
        self.instructWindow.insert('end',
            ("In the DATABASE BROWSER pane, double-click on an experiment's title; that" +
            " experiment's figure directory will immediately be added to the list.\n"),
        "level3_text")
        self.instructWindow.insert('end', "b) ", "level2_listhead")
        self.instructWindow.insert('end',
            ("In the DATABASE BROWSER pane, single-click on an experiment's title; that" +
            " experiment's figure directory will appear in the manual text entry in the" +
            " EXPERIMENT SELECTOR pane. To add it to the list, you must then confirm your" +
            " choice by pressing the \"Add\" button below the text entry.\n"),
        "level3_text")
        self.instructWindow.insert('end', "c) ", "level2_listhead")
        self.instructWindow.insert('end',
            ("In the EXPERIMENT SELECTOR pane, click the \"Browse\" button. You can then" +
            " browse the directory structure freely. When you find the figure directory" +
            " you want, single-click it and press \"Ok\"; the directory will immediately be" +
            " added to the list.\n"),
        "level3_text")
        self.instructWindow.insert('end', "d) ", "level2_listhead")
        self.instructWindow.insert('end',
            ("Manually enter the figure directory into the EXPERIMENT SELECTOR's text entry;" +
            " you must then confirm your choice by pressing the \"Add\" button.\n\n"),
        "level3_text")

        self.instructWindow.insert('end',
            ("At any time, you may remove an experiment from the EXPERIMENT SELECTOR's list by" +
            " selecting it and pressing the \"Delete\" button. You can also reorder the experiments" +
            " in the list by dragging and dropping.\n\n"),
        "level2_text")

        self.instructWindow.insert('end', "2) ", "level1_listhead")
        self.instructWindow.insert('end',
            ("Once you have selected the experiments you wish to view, click the \"List Figures\"" +
            " button in the EXPERIMENT SELECTOR pane. The FIGURE SELECTOR pane will then be" +
            " populated with all valid figures that were found in the chosen directories. Each" +
            " entry is followed by a number that indicates how many times this figure was found;" +
            " for example, if you selected 2 experiments, and they both contain a figure with a" +
            " given name, then that figure will be followed by a (2).\n\n"),
        "level2_text")

        self.instructWindow.insert('end', "3) ", "level1_listhead")
        self.instructWindow.insert('end',
            ("Once the list has been populated, simply select the figure you wish to view, and" +
            " press the \"Open in GhostView\" button. (Alternatively, you can simply double-click" +
            " the entry in the FIGURE SELECTOR to launch GhostView automatically.) If the figure" +
            " you selected appears in more than one experiment, a separate GhostView window will" +
            " be opened for each experiment.\n\n" +
            "To find a figure in the list more quickly, you can use the GoTo field to jump to" +
            " the first experiment that begins with your search term; or, you can use the Filter" +
            " field to narrow down the list to just those entries that match your search term.\n\n"),
        "level2_text")

        self.instructWindow.insert('end',
            ("Note that you can configure GhostView in advance by opening the Options menu and" +
            " selecting \"GhostView\". Here you can select the application you wish to use (either" +
            " \"gv\" or \"mgv\"), and configure its options, such as orientation, zoom level, etc." +
            " These settings go into effect immediately, and will be applied to all subsequent" +
            " GhostView windows (the changes don't affect GhostView windows that are already open).\n\n"),
        "level1_text")

        self.instructWindow.insert('end',
            ("MISCELLANEOUS NOTES\n" +
             "---------------------------------------\n\n"),
        "level1_text")

        self.instructWindow.insert('end', "* ", "level1_listhead")
        self.instructWindow.insert('end',
            ("The Figure Browser supports .PDF files as well; if you select a .PDF file," +
            " it will open in \"acroread\" with no input parameters (regardless of what options" +
            " were specified in the GhostView Options dialog).\n\n"),
        "level2_text")

        self.instructWindow.insert('end', "* ", "level1_listhead")
        self.instructWindow.insert('end',
            ("You can use command-line arguments to supply a list of experiments, which will be" +
            " pre-loaded into the Selected Experiments list on startup. This is done using the" +
            " flag -d (or --dirs), followed by a comma-delineated list of directories.\n\n"),
        "level2_text")

        self.instructWindow.insert('end', "* ", "level1_listhead")
        self.instructWindow.insert('end',
            ("Several settings can be configured by using the Options menu and selecting. Alternatively," +
            " you can \"save\" your settings by creating a small configuration file, which will be read" +
            " the next time the application is launched (This file is created for you automatically the" +
            " first time you run the application). The file must be called \".fbconfig\" and must be" +
            " placed in your home directory. It uses a 'key = value' format. Currently, the following" +
            " options are accepted:\n\n"),
        "level2_text")

        self.instructWindow.insert('end',
            ("GENERAL:\n" +
            "fontfamily = <name of font family>\n" +
            "fontsize = <integer>\n" +
            "browsedir = [<initial browse directory> | pwd]\n" +
            "archivewarning = [t | f]\n\n" +
            "GV / MGV:\n" +
            "app = [gv | mgv]\n" +
            "gridfigures = [t | f]\n" +
            "forceorder = [t | f]\n" +
            "orient = [portrait | landscape | upsidedown | seascape]\n" +
            "antialias = [t | f]\n\n" +
            "GV ONLY:\n" +
            "scale = <integer between -5 and 5, inclusive>\n\n" +
            "MGV ONLY:\n" +
            "autoheight = [t | f]\n" +
            "autowidth = [t | f]\n" +
            "magstep = <integer between -5 and 5, inclusive>\n\n"),
        "level3_text")

        self.instructWindow.insert('end',
            ("Not all options need be present in the file; any missing options will be set to" +
            " a default value by the program. Also, the rows may be in any order.\n\n"),
        "level2_text")

        self.instructWindow.insert('end', "Option Details:\n", "level2_text")

        self.instructWindow.insert('end', "- ", "level2_listhead")
        self.instructWindow.insert('end',
            "'fontfamily' is the name of a valid font family (i.e. helvetica, courier, etc.)\n",
        "level3_text")

        self.instructWindow.insert('end', "- ", "level2_listhead")
        self.instructWindow.insert('end',
            "'fontsize' is an integer representing the font size to use.\n",
        "level3_text")

        self.instructWindow.insert('end', "- ", "level2_listhead")
        self.instructWindow.insert('end',
            ("'browsedir' is the initial browse directory used by the EXPERIMENT SELECTOR's Browser." +
            " Alternatively, this value can be 'pwd', which causes the Browser to start in your" +
            " current working directory.\n"),
        "level3_text")

        self.instructWindow.insert('end', "- ", "level2_listhead")
        self.instructWindow.insert('end',
            ("'archivewarning' indicates whether you wish to be notified when you attempt to view" +
            " a directory in /archive.\n"),
        "level3_text")

        self.instructWindow.insert('end', "- ", "level2_listhead")
        self.instructWindow.insert('end',
            ("'app' is the GhostView application you wish to use. It can be either 'gv' or 'mgv'.\n"),
        "level3_text")

        self.instructWindow.insert('end', "- ", "level2_listhead")
        self.instructWindow.insert('end',
            ("'gridfigures' indicates whether you want all GhostView windows to be arranged in an" +
             " evenly-spaced grid layout on the screen. Up to 8 figures at a time can be arranged" +
             " in this manner.\n"),
        "level3_text")

        self.instructWindow.insert('end', "- ", "level2_listhead")
        self.instructWindow.insert('end',
            ("'forceorder' indicates whether you want to force all GhostView windows to open in" +
             " the order they were selected; i.e. wait for each window to appear before opening" +
             " the next one in the list. 'forceorder' can be either 't' or 'f'.\nNOTE: If" +
             " any of the selected figures are in the /archive directories, they may take a long" +
             " time to appear; hence, you may not wish to wait for these figures to appear before" +
             " opening subsequent figures. In this case, 'forceorder' should be set to 'f' (false).\n"),
        "level3_text")

        self.instructWindow.insert('end', "- ", "level2_listhead")
        self.instructWindow.insert('end',
            ("All other options are specific to GV, MGV, or both.\n\n"),
        "level3_text")

        self.instructWindow.insert('end',
            "Example .fbconfig file:\n",
        "level2_text")
        self.instructWindow.insert('end',
            ("fontfamily = helvetica\n" +
            "fontsize = 14\n" +
            "browsedir = /home/myname/\n" +
            "app = mgv\n" +
            "orient = portrait\n" +
            "antialias = t\n" +
            "magstep = 2\n" +
            "autoheight = f\n" +
            "autowidth = t\n"),
        "level3_text")

        self.instructWindow.configure(text_state = 'disabled')


    def ShowInstructions(self, *args):
        self.instructWindow.show()


    def SetupAboutDialog(self, *args):
        Pmw.aboutcontact("Developed by:\tJared Cohen\n"
                       + "With the help of:\tBruce Wyman\n\n\n"
                       + "For more information, contact:\n\n"
                       + "Jared Cohen\n"
                       + "The Geophysical Fluid Dynamics Lab\n"
                       + "Phone: (609) 452-6500 ext. 6972\n"
                       + "Email: Jared.Cohen@noaa.gov")
        self.aboutDialog = Pmw.AboutDialog(self.root,
                                           title = "About",
                                           applicationname = "Figure Browser")
        self.aboutDialog.transient(self.root)
        self.aboutDialog.withdraw()

    def ShowAboutDialog(self, *args):
        self.aboutDialog.show()

    def SetupChangeLog(self, *args):
        self.changeWindow = Pmw.TextDialog(self.root,
                                           title = "Change Log",
                                           text_wrap='word')
        self.changeWindow.withdraw()
        self.changeWindow.transient(self.root)

        self.changeWindow.insert('end',
            ("4/4/05\n" +
             "Bug Fix: Double-clicking an entry in the Database Browser worked correctly, but" +
             " selecting an entry and pressing the \"Add\" button did not. This has been fixed.\n\n"))

        self.changeWindow.insert('end',
            ("3/14/05\n" +
             "Updated layout: the Experiment Selector pane is now divided into two" +
             " subpanes; this allows the listbox and entry field to be independently" +
             " resized.\n\n"))

        self.changeWindow.insert('end',
            ("3/10/05\n" +
             "Bug Fix: The \"Save Options to .fbconfig\" command now works correctly.\n\n"))

        self.changeWindow.insert('end',
            ("3/7/05\n" +
             "-The Experiment Selector's list (containing selected experiments) now allows" +
             " drag-and-drop reordering of rows.\n" +
             "-Added 'archivewarning' option to the .fbconfig file.\n" +
             "-Added ability to save menu-selected options back into the .fbconfig file.\n" +
             "-Added command-line argument to specify experiments to be pre-loaded" +
             " into the Selected Experiments list." +
             "\n\n"))

        self.changeWindow.insert('end',
            ("1/26/05\n" +
             "Improved grid layout; program now attempts to place figures side-by-side if they can" +
             " all fit on screen without overlap.\n\n"))

        self.changeWindow.insert('end',
            ("1/24/05\n" +
             "Bug Fix: The grid layout is now an option; it can be set in the .fbconfig file using" +
             " the 'gridfigures' option.\n\n"))

        self.changeWindow.insert('end',
            ("1/24/05\n" +
             "Bug Fix: The new grid layout was preventing mgv from displaying multiple windows." +
             " This has been corrected.\n\n"))

        self.changeWindow.insert('end',
            ("1/19/05\n" +
             "Added an option to force GhostView windows to appear onscreen in the order" +
             " they were selected. This corresponds to the 'forceorder' option in the" +
             " .fbconfig file.\n\n"))

        self.changeWindow.insert('end',
            ("1/19/05\n" +
             "GhostView windows are now arranged in an evenly-spaced grid layout on the screen.\n\n"))

        self.changeWindow.insert('end',
            ("1/14/05\n" +
             "The .fbconfig file now supports GhostView options.\n\n"))

        self.changeWindow.insert('end',
            ("1/14/05\n" +
             "The .fbconfig file now uses a 'key = value' format to define its options." +
             " The options can now be in any order, and the file does not need to contain" +
             " all possible options (some may be omitted).\n\n"))

        self.changeWindow.insert('end',
            ("1/12/05\n" +
             "The program now recognizes 'pwd' as a valid 3rd line in the .fbconfig file." +
             " Using this value will set the initial browsing directory to your current" +
             " working directory.\n\n"))

        self.changeWindow.insert('end',
            ("1/10/05\n" +
             "The Database Browser and Figure Selector panes now" +
             " use a collapsible tree structure.\n\n"))

        self.changeWindow.insert('end',
            ("1/4/05\n" +
             "Initial release."))

        self.changeWindow.configure(text_state = 'disabled')

    def ShowChangeLog(self, *args):
        self.changeWindow.show()


    def SetupArchiveWarning(self, *args):
        self.archiveWarningWindow = Tkinter.Toplevel()
        self.archiveWarningWindow.withdraw()
        self.archiveWarningWindow.title("WARNING")
        self.archiveWarningWindow.transient(self.root)

        warningText = Tkinter.Text(self.archiveWarningWindow,
                      wrap='word',
                      width=40,
                      height=10)
        warningText.tag_config("all", justify='center')
        warningText.insert('end',
                            ("One or more of the selected directories are" +
                             " located in /archive. Figures in these" +
                             " directories may take a substantial amount of" +
                             " time to open in GhostView. Please ensure that" +
                             " all archived files have been migrated to disk" +
                             " before attempting to open them in GhostView."),
                           "all")
        warningText.grid(row=0, rowspan=3, column=0, columnspan=3, sticky='nsew')
        warningText.config(state='disabled')

        makeResizable(self.archiveWarningWindow)


        def ToggleArchiveWarningVar(*args):
            self.archiveWarningVar.set(not self.dontShowVar.get())
        def ToggleDontShowVar(*args):
            self.dontShowVar.set(not self.archiveWarningVar.get())

        self.dontShowVar = Tkinter.BooleanVar()
        self.dontShowVar.set(not self.archiveWarningVar.get())
        self.dontShowVar.trace_variable('w', ToggleArchiveWarningVar)
        self.archiveWarningVar.trace_variable('w', ToggleDontShowVar)


        dontShowCheckbutton = Tkinter.Checkbutton(self.archiveWarningWindow,
                                                  text="Don't show this warning again",
                                                  font=self.entryFieldFont,
                                                  variable=self.dontShowVar)
        dontShowCheckbutton.grid(row=3, rowspan=1, column=0, columnspan=3, sticky='ew')

        okButton = Tkinter.Button(self.archiveWarningWindow,
                                  text = "Ok",
                                  command=self.archiveWarningWindow.withdraw)
        okButton.grid(row=4, rowspan=1, column=1, columnspan=1, sticky='nsew')

        self.root.update()
        width = self.archiveWarningWindow.winfo_reqwidth()
        height = self.archiveWarningWindow.winfo_reqheight()
        self.archiveWarningWindow.geometry(str(width) + "x"
                                         + str(height) + "+"
                                         + str((self.root.winfo_screenwidth() - width)/2) + "+"
                                         + str((self.root.winfo_screenheight() - height)/2))



    def ShowFigDir(self, *args):
        if self.dbFrameSelected.get() == False:
            return
        self.dbFrameSelected.set(False)

        dir = self.dbFrame.dir
        if (dir == "") or (dir == "DIRECTORY"):
            return
        self.expBrowser.addExp_Entry.setvalue(dir)

    def AddFigDir(self, *args):
        if self.dbFrameConfirmed.get() == False:
            return
        self.dbFrameConfirmed.set(False)

        dir = self.dbFrame.dir
        if dir == "":
            self.root.bell()
            self.messageBar.message("userevent",
                                    "The selected experiment has no figures. Please select again.")
            return
        elif dir == "DIRECTORY":
            #self.messageBar.message("userevent",
            #                        ("This is a directory of experiments; you must select"
            #                        + " an individual experiment."))
            return
        self.expBrowser.AddExperiment(dir)

    def ShowFileTree(self, *args):
        if self.expBrowserModified.get() == False:
            return
        self.expBrowserModified.set(False)

        dirs = self.expBrowser.selExps
        self.figBrowser.ShowFileTree(dirs)

        if self.archiveWarningVar.get():
            inArchive = False
            for d in dirs:
                if d.startswith("/archive/"):
                    inArchive = True
                    break
            if inArchive:
                self.archiveWarningWindow.deiconify()

    def LaunchGhostView(self, *args):
        if self.figBrowserModified.get() == False:
            return
        self.figBrowserModified.set(False)

        figures = self.figBrowser.selFigs

        if figures[0].lower().endswith(".pdf"):
            self.messageBar.message("userevent",
                                    "Attempting to launch Acrobat Reader; this may take a few moments.")
            app = "acroread"
        else:
            self.messageBar.message("userevent",
                                    "Attempting to launch GhostView; this may take a few moments.")
            app = self.appVar.get()

        self.root.update_idletasks()

        command = app

        if app != "acroread":
            for key, var in self.opt_controlVars.iteritems():
                command = command + " " + self.opt_prefixes[key] + str(var.get())

        command = command + " "
        numFigs = len(figures)

        if (numFigs <= 1 or numFigs > 8) or (self.gridFiguresVar.get() == False):
            for f in range(len(figures)):
                fig = figures[f]
                os.system(command + fig + " &")

                #wait for this figure to appear before opening the next one
                if self.forceOrderVar.get():
                    lines = []
                    while lines == []:
                        (stdout, stdin, stderr) = os.popen3("ps -C " + app)
                        lines = stdin.readlines()[1:]
                    pid = (lines[-1].split(" "))[0]

                    lines = []
                    while lines == []:
                        if app == "gv":
                            (stdout, stdin, stderr) = os.popen3("xwininfo -name \""
                                                                + app + ": " + fig + "\"")
                        elif app == "mgv":
                            (stdout, stdin, stderr) = os.popen3("xwininfo -name \""
                                                                + "MGv(" + pid + "-1):  " + fig + "\"")
                        else: #app == "acroread"
                            fileName = fig.split("/")[-1]
                            (stdout, stdin, stderr) = os.popen3("xwininfo -name \"" + fileName + "\"")

                        lines = stdin.readlines()

                    stdout.close()
                    stdin.close()
                    stderr.close()
        else:
            #display the first window
            fig = figures[0]
            if app == "gv":
                os.system(command + fig + " -geometry +0+0 -nodsc" + " &")
            elif app == "mgv":
                os.system(command + fig + " -geometry +0+0" + " &")
            else: #app == "acroread"
                os.system(command + "-geometry +0+0 " + fig + " &")

#-->cjg: This piece of code is causing fbrowser to freeze on RHEL 5. This code
#        is related to the option "Display figures in grid layout on screen".
#        Commented out for now.
#
#            #get its geometry
#            lines = []
#            while lines == []:
#                (stdout, stdin, stderr) = os.popen3("ps -C " + app)
#                lines = stdin.readlines()[1:]
#
#            pid = (lines[-1].split(" "))[0]
#
#            lines = []
#            while lines == []:
#                if app == "gv":
#                    (stdout, stdin, stderr) = os.popen3("xwininfo -name \""
#                                                        + app + ": " + fig + "\"")
#                elif app == "mgv":
#                    (stdout, stdin, stderr) = os.popen3("xwininfo -name \""
#                                                        + "MGv(" + pid + "-1):  " + fig + "\"")
#                else: #app == "acroread"
#                    fileName = fig.split("/")[-1]
#                    (stdout, stdin, stderr) = os.popen3("xwininfo -name \"" + fileName + "\"")
#
#                lines = stdin.readlines()
#
#            stdout.close()
#            stdin.close()
#            stderr.close()
#
#            width_line = lines[7].strip(" \n")
#            height_line = lines[8].strip(" \n")
#            width = int(width_line.split(" ")[1])
#            height = int(height_line.split(" ")[1])
#
#            #find the left border of the rightmost figure
#            availWidth = self.scr_width - width
#
#            #calculate the layout of the other figures
#
#            topHalf = int(math.ceil(numFigs/2.0)) #the number of figures in the top half of the screen
#                                                  # (if any)
#
#            fullHeight = (height >= self.scr_height)  #does one figure take up the screen's full height?
#            fitVert = (2 * height <= self.scr_height)    #will all figures fit vertically on screen?
#            fitHoriz = (topHalf * width <= self.scr_width) #will all figures fit horizontally on screen?
#
#            right = str(self.scr_width - width)
#            bottom = str(self.scr_height - height)
#
#            geoms = []
#            if numFigs == 2:
#                geoms.append("+0+0")
#                if fitHoriz:
#                    geoms.append("+" + str(width) + "+0")
#                else:
#                    geoms.append("+" + right + "+0")
#            else:
#                if fitHoriz and fitVert:
#                    for i in range(0, topHalf):
#                        geoms.append("+" + str(i*width) + "+0")
#                    for i in range(0, (numFigs-topHalf)):
#                        geoms.append("+" + str(i*width) + "+" + str(height))
#
#                elif fitVert:
#                    for i in range(0, topHalf):
#                        geoms.append("+" + str((i*availWidth) / (topHalf-1)) + "+0")
#                    for i in range(0, (numFigs-topHalf)):
#                        geoms.append("+" + str((i*availWidth) / (topHalf-1)) + "+" + str(height))
#
#                elif fitHoriz:
#                    if not fullHeight:
#                        for i in range(0, topHalf):
#                            geoms.append("+" + str(i*width) + "+0")
#                        for i in range(0, (numFigs-topHalf)):
#                            geoms.append("+" + str(i*width) + "+" + bottom)
#
#                    else:
#                        for i in range(numFigs):
#                            geoms.append("+" + str(i*width) + "+0")
#
#                else:
#                    if not fullHeight:
#                        for i in range(0, topHalf):
#                            geoms.append("+" + str((i*availWidth) / (topHalf-1)) + "+0")
#                        for i in range(0, (numFigs-topHalf)):
#                            geoms.append("+" + str((i*availWidth) / (topHalf-1)) + "+" + bottom)
#
#                    else:
#                        for i in range(numFigs):
#                            geoms.append("+" + str((i*availWidth)/(numFigs-1)) + "+0")
#<--cjg

            #display the other figures
            for f in range(1, len(figures)):
                fig = figures[f]

                if app == "gv":
# cjg               os.system(command + fig + " -geometry " + geoms[f] + " -nodsc &")
                    os.system(command + fig + " -nodsc &")
                elif app == "mgv":
# cjg               os.system(command + fig + " -geometry " + geoms[f] + " &")
                    os.system(command + fig + " &")
                else: #app == "acroread"
# cjg               os.system(command + "-geometry " + geoms[f] + " " + fig + " &")
                    os.system(command + "-geometry " + fig + " &")

#-->cjg: This piece of code is causing fbrowser to freeze on RHEL 5. This code
#        is related to the option "Force figures to appear in order of selection".
#        Commented out for now.
#               #wait for this figure to appear before opening the next one
#               if self.forceOrderVar.get():
#                   lines = []
#                   while lines == []:
#                       (stdout, stdin, stderr) = os.popen3("ps -C " + app)
#                       lines = stdin.readlines()[1:]
#                   pid = (lines[-1].split(" "))[0]

#                   lines = []
#                   while lines == []:
#                       if app == "gv":
#                           (stdout, stdin, stderr) = os.popen3("xwininfo -name \""
#                                                               + app + ": " + fig + "\"")
#                       elif app == "mgv":
#                           (stdout, stdin, stderr) = os.popen3("xwininfo -name \""
#                                                               + "MGv(" + pid + "-1):  " + fig + "\"")
#                       else: #app == "acroread"
#                           fileName = fig.split("/")[-1]
#                           (stdout, stdin, stderr) = os.popen3("xwininfo -name \"" + fileName + "\"")

#                       lines = stdin.readlines()

#                   stdout.close()
#                   stdin.close()
#                   stderr.close()
#<--cjg


########################################################################

def main():

    root = Tix.Tk()
    Pmw.initialise(root)
    browser = fbMainGui(root)
    root.mainloop()

if __name__=='__main__': main()
