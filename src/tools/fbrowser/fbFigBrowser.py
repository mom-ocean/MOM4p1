#!/usr/bin/env python

import sys, Tkinter, tkFileDialog, types, tkFont, Pmw, Tix
import distutils.filelist as filelist

__FILE__ = sys._getframe().f_code.co_filename

def makeResizable(widget):
    for i in range(0, widget.grid_size()[0]):
        widget.columnconfigure(i, weight=1)
    for i in range(0, widget.grid_size()[1]):
        widget.rowconfigure(i, weight=1)

extensions = [".ps", ".ps.z", ".ps.gz", ".ps.old", ".pdf"]


class MyTreeNode:
    def __init__(self, parent, key, fullpath="", show=True):
        self.parent = parent
        self.key = key
        self.fullpaths=[fullpath,]
        self.defaultShow=show
        self.show=show
        self.index=-1

        if self.parent and self.parent.relPath != "":
            if not self.parent.relPath.endswith("/"):
                self.relPath = self.parent.relPath + "/" + self.key
            else:
                self.relPath = self.parent.relPath + self.key
        else:
            self.relPath = self.key

        self.children = []
        self.refCount = 0

        #print "New TreeNode: key = %s, fullpath = %s, relpath = %s" % (self.key, fullpath, self.relPath)

    def GetChild(self, key):
        for c in self.children:
            if c.key == key:
                return c
        return None

    def AddChild(self, key, fullPath, show=True):
        child = MyTreeNode(parent=self, key=key, fullpath = fullPath, show=show)
        self.children.append(child)
        return child

    def SortChildren(self):
        for i in range(0, len(self.children)-1):
            min = i
            for j in range((i+1), len(self.children)):
                if self.children[j].key < self.children[min].key:
                    min = j
            temp = self.children[i]
            self.children[i] = self.children[min]
            self.children[min] = temp

        files = []
        subdirs = []
        for c in self.children:
            if len(c.children) > 0:
                subdirs.append(c)
            else:
                files.append(c)
            c.refCount = len(c.fullpaths)

        self.children = subdirs
        self.children.extend(files)


class MyTree:
    def __init__(self, key, showRoot=True):
        self.rootNode = MyTreeNode(parent=None, key=key, fullpath=key, show=showRoot)

    def AddNode(self, relPath, fullPath, show=True):
        curNode = self.rootNode

        relPathStr = "/".join(relPath)
        prevPath = fullPath[:(len(fullPath)-len(relPathStr))]

        for i in range(len(relPath)):
            curKey = relPath[i]
            node = curNode.GetChild(curKey)

            if curNode != self.rootNode and not (prevPath in curNode.fullpaths):
                curNode.fullpaths.append(prevPath)

            if node != None:
                curNode = node
            elif curKey != None and curKey != "":
                if len(prevPath) > 0 and prevPath[-1] == "/":
                    node = curNode.AddChild(curKey, (prevPath + curKey), show)
                else:
                    node = curNode.AddChild(curKey, (prevPath + "/" + curKey), show)
                curNode = node

            if len(prevPath) > 0 and prevPath[-1] == "/":
                prevPath = prevPath + curKey
            else:
                prevPath = prevPath + "/" + curKey

        if curNode != self.rootNode and not (prevPath in curNode.fullpaths):
            curNode.fullpaths.append(prevPath)

    def RemoveBadFiles(self, curNode):
        if len(curNode.children) == 0: #the node is a file
            for ext in extensions:
                if curNode.key.lower().endswith(ext):
                    return
            curNode.parent.children.remove(curNode)
            del curNode
            return
        else: #the node is a directory
            i = len(curNode.children)-1
            while i >= 0:
                self.RemoveBadFiles(curNode.children[i])
                i -= 1


class fbFigBrowser:

    def __init__(self, root=None, modified=None, messageBar=None,
                 gvButtonFont=None, entryFieldFont=None):
        """
        Constructor
        root: the parent container widget (e.g. Toplevel, Frame, etc.)
        modified: a Tkinter.BooleanVar control variable
        messageBar: a Pmw.MessageBar object
        gvButtonFont: a tkFont.Font object
        entryFont: a tkFont.Font object
        """
        if root:
            self.root = root
        else:
            self.root = Tix.Tk()

        self.messageBar = messageBar

        self.modified = modified
        if self.modified is None:
            self.modified = Tkinter.BooleanVar()
            self.modified.set(False)

        self.gvButtonFont = gvButtonFont
        if self.gvButtonFont is None:
            self.gvButtonFont = tkFont.Font(family="Arial", size=14, weight="bold")

        self.entryFieldFont = entryFieldFont
        if self.entryFieldFont is None:
            self.entryFieldFont = tkFont.Font(family="helvetica", size=12, weight="normal")

        self.frame = Tkinter.Frame(self.root)

        self.fileList = []
        self.nodeList = []
        self.indexList = []
        self.expandList = []

        self.goto_EntryField = Pmw.EntryField(self.frame,
                                              modifiedcommand=self.FindFirstEntry,
                                              entry_font=self.entryFieldFont,
                                              entry_bg="#d8d8d8",
                                              label_text="Go To:",
                                              labelpos="w")
        self.goto_EntryField.grid(row=0, rowspan=1, column=0, columnspan=5, sticky='sew')

        self.filter_EntryField = Pmw.EntryField(self.frame,
                                                command=self.FilterEntries,
                                                entry_font=self.entryFieldFont,
                                                entry_bg="#d8d8d8",
                                                label_text="Filter:",
                                                labelpos="w")
        self.filter_EntryField.grid(row=0, rowspan=1, column=6, columnspan=5, sticky='sew')

        self.filter_Button = Tkinter.Button(self.frame,
                                            command=self.FilterEntries,
                                            text="Search")
        self.filter_Button.grid(row=0, rowspan=1, column=11, columnspan=1, sticky='sew')

        self.figs_Tree = Tix.Tree(self.frame, command=self.GhostView)
        self.figs_Tree.hlist.configure(separator = "/",
                                       selectmode = "single")
        self.figs_Tree.grid(row=1, rowspan=9, column=0, columnspan=12, sticky='nsew')


        self.expandAll_Button = Tkinter.Button(self.frame,
                                               command=self.ExpandAll,
                                               text="Expand All")
        self.expandAll_Button.grid(row=10, rowspan=1, column=0, columnspan=1, sticky='new')

        self.collapseAll_Button = Tkinter.Button(self.frame,
                                                 command=self.CollapseAll,
                                                 text="Collapse All")
        self.collapseAll_Button.grid(row=10, rowspan=1, column=1, columnspan=1, sticky='new')

        self.openButton = Tkinter.Button(self.frame,
                                         text="Open in GhostView",
                                         font = self.gvButtonFont,
                                         command=self.GhostView)
        self.openButton.grid(row=10, rowspan=1, column=10, columnspan=2, sticky='nsew')

        self.root.update()


    def ShowFileTree(self, dirs):
        self.BuildFileTree(dirs)
        self.figs_Tree.hlist.delete_all()
        self.BuildTreeWidget(self.tree.rootNode)

    def BuildFileTree(self, dirs):
        self.tree = MyTree("", showRoot=False)
        files = []
        for i in range(len(dirs)):
            dir = dirs[i]
            self.tree.AddNode([dir,], dir, False)
            files.append(list(filelist.findall(dir)))

            relFiles = []
            for j in range(len(files[i])):
                f = (files[i])[j]
                relFiles.append(list(f[len(dir):].strip("/").split("/")))
                self.tree.AddNode(relFiles[j], (files[i])[j], True)

        self.tree.RemoveBadFiles(self.tree.rootNode)

    def BuildTreeWidget(self, curNode):
        self.fileList = []
        self.nodeList = []
        self.Preorder_BuildTreeWidget(curNode)
        self.figs_Tree.autosetmode()

        self.indexList = range(len(self.nodeList))
        for i in range(len(self.nodeList)):
            self.nodeList[i].index = i
            self.expandList.append(False)
            self.figs_Tree.close(self.nodeList[i].relPath)

    def Preorder_BuildTreeWidget(self, curNode):
        if curNode.show:
            self.figs_Tree.hlist.add(curNode.relPath,
                                     text = ((curNode.key + "  (%d)") % curNode.refCount) )
            self.fileList.append(list(curNode.fullpaths))
            self.nodeList.append(curNode)
        curNode.SortChildren()
        for c in curNode.children:
            self.Preorder_BuildTreeWidget(c)

    def RebuildTreeWidget(self, curNode):
        self.figs_Tree.hlist.delete_all()
        self.indexList = []
        self.Preorder_RebuildTreeWidget(curNode)
        self.figs_Tree.autosetmode()
        for i in range(len(self.nodeList)):
            if self.nodeList[i].show:
                if self.expandList[i]:
                    self.figs_Tree.open(self.nodeList[i].relPath)
                else:
                    self.figs_Tree.close(self.nodeList[i].relPath)

    def Preorder_RebuildTreeWidget(self, curNode):
        if curNode.show:
            self.indexList.append(curNode.index)
            self.figs_Tree.hlist.add(curNode.relPath,
                                     text = ((curNode.key + "  (%d)") % curNode.refCount) )
        for c in curNode.children:
            self.Preorder_RebuildTreeWidget(c)

    def FindFirstEntry(self, *args):
        if len(self.indexList) == 0:
            self.root.bell()
            if self.messageBar:
                self.messageBar.message("userevent", "The figure list is empty.")
            return
        text = self.goto_EntryField.getvalue().lower()
        if text == "":
            self.figs_Tree.hlist.yview(self.nodeList[self.indexList[0]].relPath)
            return
        goto_index = -1
        for i in range(len(self.indexList)):
            index = self.indexList[i]
            if self.nodeList[index].key.lower().startswith(text) and\
                            self.nodeList[index].parent == self.tree.rootNode:
                goto_index = i
                break
        if goto_index == -1:
            return

        relPath = self.nodeList[self.indexList[goto_index]].relPath
        self.root.update()

        self.figs_Tree.hlist.yview(relPath)

    def FilterEntries(self, *args):
        if len(self.nodeList) == 0:
            self.root.bell()
            if self.messageBar:
                self.messageBar.message("userevent", "The figure list is empty.")
            return
        text = self.filter_EntryField.getvalue().lower()

        for i in range(len(self.nodeList)):
            self.nodeList[i].show = self.nodeList[i].defaultShow
            self.expandList[i] = False
        if text != "":
            self.rec_FilterEntries(self.tree.rootNode, text, False)

        self.RebuildTreeWidget(self.tree.rootNode)

        if self.goto_EntryField.getvalue() != "":
            self.FindFirstEntry()

    def rec_FilterEntries(self, curNode, text, foundParent):
        found_me = (curNode.key.lower().find(text) != -1)
        found_parent = foundParent
        found_child = False

        for c in curNode.children:
            if self.rec_FilterEntries(c, text, found_me):
                found_child = True

        curNode.show = curNode.defaultShow and (found_me or found_parent or found_child)
        self.expandList[curNode.index] = found_child

        return (found_me or found_child)

    def ExpandAll(self, *args):
        if len(self.nodeList) == 0:
            self.root.bell()
            if self.messageBar:
                self.messageBar.message("userevent", "The figure list is empty.")
            return
        for i in range(len(self.indexList)):
            path = self.nodeList[self.indexList[i]].relPath
            self.figs_Tree.open(path)

    def CollapseAll(self, *args):
        if len(self.nodeList) == 0:
            self.root.bell()
            if self.messageBar:
                self.messageBar.message("userevent", "The figure list is empty.")
            return
        for i in range(len(self.indexList)):
            path = self.nodeList[self.indexList[i]].relPath
            self.figs_Tree.close(path)





    def GhostView(self, *args):
        curSel = self.figs_Tree.hlist.info_selection()
        if len(curSel) == 0:
            self.root.bell()
            if self.messageBar:
                self.messageBar.message("userevent",
                                        "You did not select a figure to be opened in GhostView.")
            return
        curSel = curSel[0]
        line = 0
        for i in range(len(self.nodeList)):
            if self.nodeList[i].relPath == curSel:
                line = i
                break

        if len(self.nodeList[line].children) != 0:
            return
        self.selFigs = self.fileList[line]
        self.modified.set(True)



    def grid(self, **options):
        self.frame.grid(**options)
        makeResizable(self.frame)
        self.frame.rowconfigure(0, weight=0)
        self.frame.rowconfigure(10, weight=0)


########################################################################

def main():
    root = Tix.Tk()
    Pmw.initialise(root)
    #root.geometry("500x500+0+0")

    browser = fbFigBrowser(root)
    browser.frame.grid(sticky='nsew')
    makeResizable(root)

    root.mainloop()

if __name__=='__main__': main()
