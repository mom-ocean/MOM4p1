#!/usr/bin/env python

import sys, Tkinter, tkFileDialog, types, tkFont, Pmw, Tix
import MySQLdb

__FILE__ = sys._getframe().f_code.co_filename

def makeResizable(widget):
    for i in range(0, widget.grid_size()[0]):
        widget.columnconfigure(i, weight=1)
    for i in range(0, widget.grid_size()[1]):
        widget.rowconfigure(i, weight=1)

def tokensort(L):

    #split into numeric and non-numeric substrings (tokens)
    tokens = []
    for s in L:
        #start a new list of tokens
        V = []
        V.append(s[0])
        isPrevNum = False
        isCurNum = s[0].isdigit()

        #iterate through the current string
        for i in range(1, len(s)):
            isPrevNum = isCurNum
            isCurNum = s[i].isdigit()
            if isCurNum == isPrevNum: #we're still in the same token
                V[-1] = V[-1] + s[i]
            else: #start a new token
                V.append(s[i])
        tokens.append(V)

    #find the maximum number of decimal places we need
    maxlen = 0
    for i in range(len(tokens)):
        for j in range(len(tokens[i])):
            if (tokens[i][j].isdigit() and len(tokens[i][j]) > maxlen):
                maxlen = len(tokens[i][j])

    #prepend zeroes to all numeric tokens
    for i in range(len(tokens)):
        for j in range(len(tokens[i])):
            if tokens[i][j].isdigit():
                tokens[i][j] = tokens[i][j].zfill(maxlen)

    #recombine the tokens, and form a list of (new string, original string) tuples
    tupleList = []
    for i in range(len(tokens)):
        key = ("".join(tokens[i])).upper()
        val = L[i]
        tupleList.append((key, val))

    #sort the list of tuples, and return the sorted original list
    tupleList.sort()
    L[:] = [val for (key, val) in tupleList]


class fbDatabaseFrame:

    def __init__(self, root=None, selected=None, confirmed=None, messageBar=None):
        """
        Constructor
        root: the parent container widget (e.g. Toplevel, Frame, etc.)
        selected: a Tkinter.BooleanVar control variable
        confirmed: a Tkinter.BooleanVar control variable
        messageBar: a Pmw.MessageBar object
        """
        if root:
            self.root = root
        else:
            self.root = Tix.Tk()

        self.messageBar = messageBar

        self.selected = selected
        if self.selected == None:
            self.selected = Tkinter.BooleanVar()
            self.selected.set(False)

        self.confirmed = confirmed
        if self.confirmed == None:
            self.confirmed = Tkinter.BooleanVar()
            self.confirmed.set(False)

        self.frame = Tkinter.Frame(self.root)

        self.db_Tree = Tix.Tree(self.frame,
                                browsecmd = self.SelectRow,
                                command = self.ConfirmRow)
        self.db_Tree.hlist.configure(separator = "/",
                                     selectmode = "single")
        self.db_Tree.grid(row=0, rowspan=40, column=0, columnspan=10, sticky='nsew')

        self.confirmButton = Tkinter.Button(self.frame,
                                            text="Add Figure Directory",
                                            command=self.ConfirmRow)
        self.confirmButton.grid(row=40, rowspan=1, column=5, columnspan=1, sticky='nsew')

        self.LoadDatabase()
        #self.db_Tree.autosetmode()

        self.root.update()


    def LoadDatabase(self):
        db = MySQLdb.connect(db="model_development", host="cobweb", user="gfdl", passwd="wrks4me")
        db.query("show tables;")
        res = db.store_result()
        rows = list(res.fetch_row(maxrows=0))
        for r in range(len(rows)):
            rows[r] = (rows[r])[0]
        tokensort(rows)

        self.figDirs = {}

        for tableName in rows:
            self.figDirs[tableName] = "DIRECTORY"
            self.db_Tree.hlist.add(tableName, text=tableName)
            #self.db_Tree.setmode(tableName, 'close')

            db.query("select exp, diag_figs from " + tableName + " order by exp;")
            res1 = db.store_result()
            rowList = list(res1.fetch_row(maxrows=0))
            expList = []
            figDict = {}
            for s in range(len(rowList)):
                if rowList[s][1].endswith("\n"):
                    figDict[rowList[s][0]] = rowList[s][1][:-1]
                else:
                    figDict[rowList[s][0]] = rowList[s][1]
                expList.append(rowList[s][0])
            tokensort(expList)

            for s in expList:
                path = tableName + "/" + s
                self.db_Tree.hlist.add(path, text=s)
                #self.db_Tree.setmode(path, 'none')
                self.figDirs[path] = figDict[s]

        self.db_Tree.autosetmode()
        for tableName in rows:
            self.db_Tree.close(tableName)

        db.close()

    def SelectRow(self, entry):
        self.dir = self.figDirs[entry]
        self.selected.set(True)

    def ConfirmRow(self, entry=None):
        if entry == None:
            entry = self.db_Tree.hlist.info_selection()[0]
        if entry == "":
            self.root.bell()
            if self.messageBar:
                self.messageBar.message("userevent",
                                        "You did not select an experiment from the list.")
            return
        self.dir = self.figDirs[entry]
        self.confirmed.set(True)

    def grid(self, **options):
        self.frame.grid(**options)
        makeResizable(self.frame)


########################################################################

def main():

    root = Tix.Tk()
    Pmw.initialise(root)
    #root.geometry("500x500+0+0")

    frame = Tkinter.Frame(root)

    panes = Pmw.PanedWidget(root,
                            #hull_borderwidth=1,
                            hull_relief='sunken',
                            orient='horizontal',
                            handlesize=0,
                            separatorrelief='raised',
                            separatorthickness=8)
    leftPane = panes.add("leftPane", size=0.25)
    rightPane = panes.add("rightPane", size=0.75)

    dbFrame = fbDatabaseFrame(leftPane)
    dbFrame.frame.grid(sticky='nsew')
    makeResizable(root)
    makeResizable(dbFrame.frame)

    root.mainloop()

if __name__=='__main__': main()
