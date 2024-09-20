Base for Docker containers needing to include `circos` capabilities.

On the Debian version of `circos` the binary is installed to `/usr/bin` and thus the directories herein are simply copied
underneath that location. This is unorthodox but simplifies other aspects as the default location for all of these things is in
the same directory as the binary.

