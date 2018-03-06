#! /bin/env python3
# ------------------
# Functions
# ------------------

def textindent(text, indents=1):
    # {{{

    if not text or not isinstance(text, str):
        return ""
    jointext = ''.join(['\n'] + ['    '] * indents)
    return jointext.join(text.split('\n'))

    # }}}

# ------------------

def textcode(codelist):
    # {{{

    codeleading = '>>> '
    codesep = "\n" + codeleading
    code = codesep.join(codelist)
    return codeleading + code

    # }}}

# ------------------

def textheader(header):
    # {{{

    return '\n'.join([header, ''.join(['-']*len(header))])

    # }}}

# ------------------

def textdescription(textdict):
    # {{{

    text = []
    for key, describe in textdict.items():
        if isinstance(describe, list):
            describetext = textindent('\n'.join(describe), indents=1)
        elif isinstance(describe, str):
            describetext = describe
        else:
            pass
        text.append(' : '.join([key, describetext]))
    return '\n'.join(text)

    # }}}

# ------------------

def helpstring(describe,
               parameterdicts=None,
               returns=None,
               examplecodelists=None,
               sep='\n\n\n'):
    # {{{

    ds = [describe]
    if parameterdicts is not None:
        ds.append(
            '\n'.join(
                [
                    textheader('Parameters'),
                    textdescription(parameterdicts)
                ]
            )
        )
    else:
        pass
    if returns is not None:
        ds.append(
            '\n'.join(
                [
                    textheader('Returns'),
                    returns
                ]
            )
        )
    else:
        pass
    if examplecodelists is not None:
        ds.append(
            '\n'.join(
                [
                    textheader('Examples'),
                    textcode(examplecodelists)
                ]
            )
        )
    else:
        pass
    return sep.join(ds)

    # }}}


# ------------------
# Classes
# ------------------
class AppendHelp:
    def __init__(self, helpstr, join):
        self.string = helpstr
        self.join = join
    def __call__(self, func):
        if func.__doc__ is None:
            func.__doc__ = ''
        else:
            pass
        func.__doc__ = self.join.join([func.__doc__, self.string])
        return func


# ------------------
# EOF
# ------------------
