"""
pyvtk.xmlwrite.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Simple class to generate a well-formed XML file.

Copyright (c) 05-07-2021,  Shawn W. Walker
"""

_DEFAUL_ENCODING = "ASCII"

class XmlWriter:
    def __init__(self, filepath, addDeclaration = True):
        self.stream = open(filepath, "wb")
        self.openTag = False
        self.current = []
        if (addDeclaration): self.addDeclaration()

    def addComment(self, sstr):
        """ Adds (open and close) a single comment contained in string sstr. 
            TODO: Add a smart check for the position of the comments in the file. For now,
                  we rely on the caller.
        """
        if self.openTag: 
            self.stream.write(b">")
            self.openTag = False
        self.stream.write(b'\n<!-- ') # new line is not strictly necessary
        self.stream.write(sstr.encode(_DEFAUL_ENCODING))
        self.stream.write(b' -->')    # new line here is not necessary?
        
    def close(self):
        assert(not self.openTag)
        self.stream.close()

    def addDeclaration(self):
        self.stream.write(b'<?xml version="1.0"?>')
    
    def openElement(self, tag):
        if self.openTag: self.stream.write(b">")
        st = "\n<%s" % tag
        self.stream.write(st.encode(_DEFAUL_ENCODING))
        self.openTag = True
        self.current.append(tag)
        return self

    def closeElement(self, tag = None):
        if tag:
            assert(self.current.pop() == tag)
            if (self.openTag):
                self.stream.write(b">")
                self.openTag = False
            st = "\n</%s>" % tag
            self.stream.write(st.encode(_DEFAUL_ENCODING))
        else:
            self.stream.write(b"/>")
            self.openTag = False
            self.current.pop()
        return self

    def addText(self, text):
        if (self.openTag):
            self.stream.write(b">\n")
            self.openTag = False
        self.stream.write(text.encode(_DEFAUL_ENCODING))
        return self

    def addAttributes(self, **kwargs):
        assert (self.openTag)
        for key in kwargs:
            st = ' %s="%s"'%(key, kwargs[key])
            self.stream.write(st.encode(_DEFAUL_ENCODING))
        return self

#
