Python2 files were copied from the original sciara-project-tools directory.

The python3-provided conversion script called 2to3-3.7 was then run on them in sub-groups or individually, as I constructed the new directory layout.


Syntax: 2to3-3.7 [--no-diffs] --nobackups -w file|dir

Most common usage:
	2to3-3.7 --nobackups -w file


Header lines were changed automatically. 
#!/usr/bin/env python2.7
#!/usr/bin/env python

All now:
#!/usr/bin/env python3

Library .py files that did not have a header still don't.
