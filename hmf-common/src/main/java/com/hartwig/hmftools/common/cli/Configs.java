package com.hartwig.hmftools.common.cli;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class Configs {

    private static final Logger LOGGER = LogManager.getLogger(Configs.class);

    private Configs() {
    }

    @NotNull
    public static <E extends Enum<E>> E defaultEnumValue(@NotNull final CommandLine cmd, @NotNull final String argument,
            @NotNull final E defaultValue) throws ParseException {
        if (cmd.hasOption(argument)) {
            final String optionValue = cmd.getOptionValue(argument);
            try {
                final E value = E.valueOf(defaultValue.getDeclaringClass(), optionValue);
                if (!value.equals(defaultValue)) {
                    LOGGER.info("Using non default value {} for parameter {}", optionValue, argument);
                }

                return value;
            } catch (IllegalArgumentException e) {
                throw new ParseException("Invalid validation stringency: " + optionValue);
            }
        }

        return defaultValue;
    }

    public static boolean defaultBooleanValue(@NotNull final CommandLine cmd, @NotNull final String opt, final boolean defaultValue) {
        if (cmd.hasOption(opt)) {
            final boolean result = Boolean.parseBoolean(cmd.getOptionValue(opt));
            if (result  != defaultValue) {
                LOGGER.info("Using non default value {} for parameter {}", result, opt);
            }
            return result;
        }

        return defaultValue;
    }

    public static double defaultDoubleValue(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue) {
        if (cmd.hasOption(opt)) {
            final double result = Double.parseDouble(cmd.getOptionValue(opt));
            if (!Doubles.equal(result, defaultValue)) {
                LOGGER.info("Using non default value {} for parameter {}", result, opt);
            }
            return result;
        }

        return defaultValue;
    }

    public static int defaultIntValue(@NotNull final CommandLine cmd, @NotNull final String opt, final int defaultValue) {
        if (cmd.hasOption(opt)) {
            final int result = Integer.parseInt(cmd.getOptionValue(opt));
            if (result != defaultValue) {
                LOGGER.info("Using non default value {} for parameter {}", result, opt);
            }
            return result;
        }

        return defaultValue;
    }

    @NotNull
    public static String defaultStringValue(@NotNull final CommandLine cmd, @NotNull final String opt, final String defaultValue) {
        if (cmd.hasOption(opt)) {
            final String result = cmd.getOptionValue(opt);
            if (!result.equals(defaultValue)) {
                LOGGER.info("Using non default value {} for parameter {}", result, opt);
            }
            return result;
        }

        return defaultValue;
    }


    public static boolean containsFlag(@NotNull final CommandLine cmd, @NotNull final String opt) {
        if (cmd.hasOption(opt)) {
            LOGGER.info("Using non default {} flag", opt);
            return true;
        }
        return false;
    }

}
