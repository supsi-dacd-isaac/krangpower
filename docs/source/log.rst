Logging
=======

krangpower mantains two different logs that can be useful for debugging purposes.

In Windows, they can be found under:

.. code::

    %APPDATA%/krangpower/commands.log
    %APPDATA%/krangpower/krangpower.log

For Linux, they are created under:

.. code::

    /var/log/krangpower/commands.log
    /var/log/krangpower/krangpower.log

Log level
.........

In order to utilize the logs, it is of uttermost importance to correctly set the module-wide log level. The level has to be set with a dedicated function
exposed at the module level:

.. autofunction:: krangpower.set_log_level

The default/starting value of the logging level, when using :code:`krangpower` out of the box, is set at :code:`logging.WARNING`
(equivalent to a numerical value of 30), but can be customized under:

.. code::

    ./config/krang_config.cfg --> misc settings --> default_logging_level

.. warning::
    The commands.log and krangpower.log files are only written **if the log level is <= logging.DEBUG** (10).

.. warning::
    The log files are quite verbose and writing them often causes a *major* execution time overhead, so
    they should be left unactive if not needed.

.. tip::
    A log level <= logging.INFO (20) does not write the logs, but displays a bit more information on the console at runtime, notably
    progress bars for the computation/execution of properties/function such as :code:`krang.graph`, etc.

commands.log
............

This file contains every command passed to the underlying OpenDSS textual interface, complete with a timestamp and the result
returned by OpenDSS. The log is mantained directly by the :code:`krangpower.enhancer` submodule.

krangpower.log
..............

This file logs various activity aspect of the Krang components and :code:`PackedOpendssElements`. In particular, you can find
here every call to the OpenDSSDirect.py attributes and reports on the choices made by the components' :code:`_setparameter` method.
