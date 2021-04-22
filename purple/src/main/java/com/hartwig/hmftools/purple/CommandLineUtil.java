package com.hartwig.hmftools.purple;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CommandLineUtil {

    private static final Logger LOGGER = LogManager.getLogger(CommandLineUtil.class);

    private CommandLineUtil() {
    }

    @NotNull
    public static String defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt, @NotNull final String defaultValue) {
        return cmd.hasOption(opt) ? cmd.getOptionValue(opt) : defaultValue;
    }

    public static double defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue) {
        if (cmd.hasOption(opt)) {
            final double result = Double.parseDouble(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

    public static int defaultIntValue(@NotNull final CommandLine cmd, @NotNull final String opt, final int defaultValue) {
        if (cmd.hasOption(opt)) {
            final int result = Integer.parseInt(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }
}
