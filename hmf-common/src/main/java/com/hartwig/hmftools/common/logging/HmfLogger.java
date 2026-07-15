package com.hartwig.hmftools.common.logging;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.spi.ExtendedLogger;
import org.apache.logging.log4j.spi.ExtendedLoggerWrapper;

// Logger with a custom DEBUG_2 level: finer than DEBUG (intLevel 500), coarser than TRACE (600), so it can be enabled to show
// extra per-read detail without pulling in third-party TRACE noise (htsjdk and friends log at TRACE). The wrapper mirrors what
// log4j2's own custom-level generator produces, exposing TARS_LOGGER.debug2(...) alongside the standard level methods.
public final class HmfLogger extends ExtendedLoggerWrapper
{
    public static final Level DEBUG_2 = Level.forName("DEBUG_2", 550);

    private static final String FQCN = HmfLogger.class.getName();

    private HmfLogger(final ExtendedLogger logger)
    {
        super(logger, logger.getName(), logger.getMessageFactory());
    }

    public static HmfLogger getLogger(final Class<?> clazz)
    {
        return new HmfLogger((ExtendedLogger) LogManager.getLogger(clazz));
    }

    public static HmfLogger getLogger(final String name)
    {
        return new HmfLogger((ExtendedLogger) LogManager.getLogger(name));
    }

    // Resolves standard and hmftools custom level names. Referencing this class loads DEBUG_2 (via the static field above) so a
    // later Level.valueOf("DEBUG_2") does not throw.
    public static Level parseLevel(final String name)
    {
        return Level.valueOf(name);
    }

    public void debug2(final String message, final Object... params)
    {
        logIfEnabled(FQCN, DEBUG_2, null, message, params);
    }
}
