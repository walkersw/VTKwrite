"""
VTKwrite.xmlwrite.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Simple class to generate a well-formed XML file.

The file is written to a temporary file next to the target and only renamed to
the target name when close() succeeds, so an error part-way through never
leaves a truncated file (or clobbers an existing one) under the target name.

Copyright (c) 05-07-2021,  Shawn W. Walker
"""

import os
import uuid

_ENCODING = "utf-8"

def escape_attribute(value):
    """ Return str(value) escaped for use inside a double-quoted XML attribute. """
    s = str(value)
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;").replace("\n", "&#10;").replace("\r", "&#13;")
             .replace("\t", "&#9;"))

def sanitize_comment(text):
    """ Return text made legal inside an XML comment: it may not contain "--"
        and may not end with "-". """
    s = str(text)
    while "--" in s:
        s = s.replace("--", "- -")
    if s.endswith("-"):
        s += " "
    return s

class XmlWriter:
    def __init__(self, filepath, addDeclaration = True):
        self.filepath = os.path.abspath(filepath)
        # If the target is a symbolic link, write (and later replace) the file it points
        # to, as opening the link for writing would, instead of replacing the link itself.
        self._target = os.path.realpath(self.filepath)
        directory, base = os.path.split(self._target)
        # Unique hidden name in the same directory, so the final rename is atomic. The
        # base name is shortened so the temporary name is never longer than needed.
        # open(..., "xb") (unlike tempfile.mkstemp) gives the usual umask permissions.
        self._tmppath = os.path.join(directory, ".%s.%s.part" % (base[:40], uuid.uuid4().hex[:8]))
        self.stream = open(self._tmppath, "xb")
        self.openTag = False
        self.current = []
        self.closed = False
        if (addDeclaration): self.addDeclaration()

    # context manager: commit on success, discard on error
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is not None:
            self.abort()
        elif not self.closed:
            self.close()
        return False

    def __del__(self):
        # A writer that is garbage collected without close() can never be completed:
        # remove its temporary file instead of leaving it behind.
        if not getattr(self, "closed", True):
            self.abort()

    def _write(self, text):
        self.stream.write(text.encode(_ENCODING))

    def addComment(self, sstr):
        """ Adds (open and close) a single comment contained in string sstr.
            Characters that are illegal in XML comments ("--") are made safe.
            TODO: Add a smart check for the position of the comments in the file. For now,
                  we rely on the caller.
        """
        if self.openTag:
            self.stream.write(b">")
            self.openTag = False
        self._write("\n<!-- " + sanitize_comment(sstr) + " -->")

    def close(self):
        """ Finish the file and move it to its final name. """
        if self.openTag or self.current:
            unclosed = ", ".join(self.current)
            self.abort()
            raise RuntimeError("XML file '%s' closed with unclosed element(s): %s"
                               % (self.filepath, unclosed))
        self.stream.close()
        try:
            os.replace(self._tmppath, self._target)
        except OSError:
            self.closed = True
            try:
                os.remove(self._tmppath)
            except OSError:
                pass
            raise
        self.closed = True

    def abort(self):
        """ Discard the partially written file. The target file is left untouched. """
        if self.closed:
            return
        try:
            self.stream.close()
        finally:
            self.closed = True
            try:
                os.remove(self._tmppath)
            except OSError:
                pass

    def addDeclaration(self):
        self.stream.write(b'<?xml version="1.0"?>')

    def openElement(self, tag):
        if self.openTag: self.stream.write(b">")
        self._write("\n<%s" % tag)
        self.openTag = True
        self.current.append(tag)
        return self

    def closeElement(self, tag = None):
        if not self.current:
            raise RuntimeError("closeElement(%r) called with no open element" % (tag,))
        if tag:
            if self.current[-1] != tag:
                raise RuntimeError("closeElement(%r) does not match the open element <%s>"
                                   % (tag, self.current[-1]))
            self.current.pop()
            if (self.openTag):
                self.stream.write(b">")
                self.openTag = False
            self._write("\n</%s>" % tag)
        else:
            if not self.openTag:
                raise RuntimeError("closeElement() without a tag requires an element "
                                   "with no content; use closeElement(%r)" % self.current[-1])
            self.stream.write(b"/>")
            self.openTag = False
            self.current.pop()
        return self

    def addText(self, text):
        if (self.openTag):
            self.stream.write(b">\n")
            self.openTag = False
        self._write(text)
        return self

    def addAttributes(self, **kwargs):
        if not self.openTag:
            raise RuntimeError("addAttributes() must directly follow openElement()")
        for key in kwargs:
            self._write(' %s="%s"' % (key, escape_attribute(kwargs[key])))
        return self

#
