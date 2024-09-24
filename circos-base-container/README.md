Base for Docker containers needing to include `circos` capabilities.

On the Debian version of `circos` the binary is installed to `/usr/bin` and thus the directories herein are simply copied
underneath that location. This is unorthodox but simplifies other aspects as the default location for all of these things is in
the same directory as the binary. The fonts and data are in tarballs as they have never really changed throughout their usage in
the pipeline so we feel it reasonable not to version them individually.

Note this build is not configured to trigger or be triggered by any other tools, given the slow rate of change it is expected that
it will be manually managed.


