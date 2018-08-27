I/O with DSS files (debug)
''''''''''''''''''''''''''

:code:`krangpower` offers compatibilty with the classic DSS format of OpenDSS, mainly for purposes of back-compatibility
and debugging.

Loading
.......

:code:`Krang` exposes a :code:`Krang.from_dss` class method. It was written with the **only** scope of allowing direct use of existing dss scripts without having to
translate them to krangpower instructions; in all other cases, writing a krangpower script is to be always preferred.

.. automethod:: krangpower.Krang.from_dss
   :noindex:

Saving
......

:code:`Krang.save_dss` saves a dss from the series of instructions that were imparted through the
:code:`Krang.command` method.

.. warning::
   Modifications operated through other means than the text interface won't be included!

.. automethod:: krangpower.Krang.save_dss
   :noindex:
