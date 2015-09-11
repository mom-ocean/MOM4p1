#!/usr/bin/env python

import os, sys, Tkinter, tkFileDialog, types, tkFont, Pmw, Tix
__FILE__ = sys._getframe().f_code.co_filename

def makeResizable(widget):
    for i in range(0, widget.grid_size()[0]):
        widget.columnconfigure(i, weight=1)
    for i in range(0, widget.grid_size()[1]):
        widget.rowconfigure(i, weight=1)

class fbExpBrowser:

    def __init__(self, root=None, modified=None, messageBar=None, entryFont=None,
                 initBrowseDir="", initExpDirs=[]):
        """
        Constructor
        root: the parent container widget (e.g. Toplevel, Frame, etc.)
        modified: a Tkinter.BooleanVar control variable
        messageBar: a Pmw.MessageBar object
        entryFont: a tkFont.Font object
        initBrowseDir: initial browsing directory (a string)
        initExpDirs: list of directories to be loaded into the Selected Experiments list
        """
        if root:
            self.root = root
        else:
            self.root = Tix.Tk()

        self.messageBar = messageBar

        self.modified = modified
        if self.modified == None:
            self.modified = Tkinter.BooleanVar()
            self.modified.set(False)

        self.initBrowseDir = str(initBrowseDir)
        self.hasBrowsed = False

        self.entryFont = entryFont
        if self.entryFont is None:
            self.entryFont = tkFont.Font(family = "helvetica",
                                         size   = 12,
                                         weight = "normal")

        self.initExpDirs = initExpDirs

        self.frame = Tkinter.Frame(self.root)

        self.LR_panes = Pmw.PanedWidget(self.frame,
                                        #hull_relief='sunken',
                                        orient='horizontal',
                                        handlesize=0,
                                        separatorrelief='raised',
                                        separatorthickness=4)
        self.leftPane = self.LR_panes.add("leftPane", size=0.5)
        self.rightPane = self.LR_panes.add("rightPane", size=0.5)

        self.LR_panes.grid(row=0, rowspan=10, column=0, columnspan=10, sticky='nsew')
        self.leftPane.grid(sticky='nsew')
        self.rightPane.grid(sticky='nsew')

        self.selExps_Label = Tkinter.Label(self.leftPane,
                                           text="Selected Experiments")
        self.selExps_Label.grid(row=0, rowspan=1, column=0, columnspan=10, sticky='sew')

        self.selExps = []
        self.selExps_ListBox = Pmw.ScrolledListBox(self.leftPane,
                                                   #label_text = "Selected experiments",
                                                   #labelpos = "n",
                                                   listbox_selectmode="single",
                                                  )
        self.selExps_ListBox.grid(row=1, rowspan=3, column=0, columnspan=10, sticky='nsew')

        self.lb = self.selExps_ListBox.component('listbox')
        self.dragIndex = 0
        self.targetIndex = 0
        self.lastTargetIndex = 0
        self.listbox_bg = self.lb.cget('bg')
        self.drag_bg = "#d8d8d8"
        self.lb.bind("<ButtonPress-1>", self.ListBox_Click)
        self.lb.bind("<B1-Motion>", self.ListBox_Drag)
        self.lb.bind("<ButtonRelease-1>", self.ListBox_Drop)

        self.deleteButton = Tkinter.Button(self.leftPane,
                                           text="Delete",
                                           command=self.DeleteExperiments)
        self.deleteButton.grid(row=1, rowspan=1, column=10, columnspan=1, sticky='nw')

        #dummyFrame = Tkinter.Frame(self.leftPane)
        #dummyFrame.grid(row=3, rowspan=1, column=0, columnspan=11, sticky='nsew')


        self.addExp_Label = Tkinter.Label(self.rightPane,
                                          text="Add new experiment")
        self.addExp_Label.grid(row=0, rowspan=1, column=1, columnspan=10, sticky='sew')
        self.addExp_Entry = Pmw.EntryField(self.rightPane,
                                           #label_text="Add new experiment",
                                           #labelpos="n",
                                           entry_bg="#d8d8d8",
                                           command = (lambda *args:
                                                      self.AddExperiment(self.addExp_Entry.getvalue())),
                                           entry_font = self.entryFont)
        #self.addExp_Entry.grid(row=0, rowspan=2, column=1, columnspan=10, sticky='new')
        self.addExp_Entry.grid(row=1, rowspan=1, column=1, columnspan=10, sticky='new')
        #print self.addExp_Entry.cget("entry_bg")

        def ExpBrowse():
            if not self.hasBrowsed:
                dir = tkFileDialog.askdirectory(title='Choose a directory',
                                                initialdir=self.initBrowseDir)
                self.hasBrowsed = True
            else:
                dir = tkFileDialog.askdirectory(title='Choose a directory')

            if not dir:
                return None
            #self.addExp_Entry.setvalue(dir)
            self.AddExperiment(dir)

        self.browseButton = Tkinter.Button(self.rightPane,
                                           text="Browse",
                                           command=ExpBrowse)
        self.browseButton.grid(row=2, rowspan=1, column=10, columnspan=1, sticky='ne')

        self.addButton = Tkinter.Button(self.rightPane,
                                        text="Add",
                                        command = (lambda *args:
                                                    self.AddExperiment(self.addExp_Entry.getvalue()))
                                       )
        self.addButton.grid(row=1, rowspan=1, column=0, columnspan=1, sticky='ne')

        dummyFrame = Tkinter.Frame(self.rightPane)
        dummyFrame.grid(row=3, rowspan=1, column=0, columnspan=11, sticky='nsew')

        self.submitButton = Tkinter.Button(self.frame,
                                           text="List Figures",
                                           command = self.SubmitExperiments)
        self.submitButton.grid(row=10, rowspan=1, column=4, columnspan=2, sticky='nsew')

        #makeResizable(self.frame)

        if len(self.initExpDirs) > 0:
            self.invalidDir = False
            for d in self.initExpDirs:
                self.AddExperiment(d)
            if self.invalidDir:
                self.root.bell()
                if self.messageBar:
                    self.messageBar.message("userevent",
                                            ("One or more of the supplied directories were invalid."))

        makeResizable(self.leftPane)
        self.leftPane.rowconfigure(0, weight=0)
        self.leftPane.rowconfigure(1, weight=0)
        self.leftPane.columnconfigure(10, weight=0)
        makeResizable(self.rightPane)
        self.rightPane.rowconfigure(0, weight=0)
        self.rightPane.rowconfigure(1, weight=0)
        self.rightPane.columnconfigure(0, weight=0)
        self.root.update()

    def ListBox_Click(self, event):
        self.dragIndex = self.lb.nearest(event.y)
        self.targetIndex = self.dragIndex
        self.lastTargetIndex = self.targetIndex

    def ListBox_Drag(self, event):
        new_targetIndex = self.lb.nearest(event.y)
        if new_targetIndex != self.targetIndex:
            self.lastTargetIndex = self.targetIndex
            self.targetIndex = new_targetIndex
            self.lb.itemconfigure(self.targetIndex, bg=self.drag_bg)
            if self.lastTargetIndex != self.dragIndex:
                self.lb.itemconfigure(self.lastTargetIndex, bg=self.listbox_bg)

    def ListBox_Drop(self, event):
        if self.targetIndex != self.dragIndex:
            self.lb.itemconfigure(self.targetIndex, bg=self.listbox_bg)
            dragLine = self.lb.get(self.dragIndex)
            self.lb.delete(self.dragIndex)
            self.lb.insert(self.targetIndex, dragLine)
            self.lb.selection_set(self.targetIndex)

    def AddExperiment(self, val):
        if val == None or val == "":
            self.root.bell()
            if self.messageBar:
                self.messageBar.message("userevent",
                                        "You did not enter a directory.")
            return
        elif not os.path.exists(val):
            self.invalidDir = True
            return
        items = list(self.selExps_ListBox.get())
        if items.count(val) == 0:
            items.append(val)
            self.selExps_ListBox.setlist(items)
            self.addExp_Entry.clear()
        else:
            self.root.bell()
            if self.messageBar:
                self.messageBar.message("userevent",
                                        "This experiment is already in the list. Please select again.")

    def DeleteExperiments(self, *args):
        selected = list(self.selExps_ListBox.curselection())
        for i in selected:
            self.selExps_ListBox.delete(i)

    def SubmitExperiments(self, *args):
        self.selExps = list(self.selExps_ListBox.get())
        if len(self.selExps) == 0:
            self.root.bell()
            if self.messageBar:
                self.messageBar.message("userevent",
                                        "No experiments to submit.")
            return
        self.modified.set(True)

    def grid(self, **options):
        self.frame.grid(**options)
        makeResizable(self.frame)
        self.frame.rowconfigure(10, weight=0)


########################################################################

def main():

    root = Tix.Tk()
    Pmw.initialise(root)
    #root.geometry("500x500+0+0")

    browser = fbExpBrowser(root)
    browser.frame.grid(sticky='nsew')
    makeResizable(root)

    root.mainloop()

if __name__=='__main__': main()
