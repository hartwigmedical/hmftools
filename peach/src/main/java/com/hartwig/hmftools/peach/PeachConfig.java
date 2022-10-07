package com.hartwig.hmftools.peach;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class PeachConfig {
    public PeachConfig(final CommandLine cmd)
    {
    }
    public static Options createOptions()
    {
        final Options options = new Options();
        return options;
    }
}
