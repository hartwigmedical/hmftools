package com.hartwig.hmftools.common.logging;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.Level;
import org.junit.Test;

public class HmfLoggerTest
{
    @Test
    public void debug2RegisteredAndParseable()
    {
        assertEquals(HmfLogger.DEBUG_2, HmfLogger.parseLevel("DEBUG_2"));
        assertEquals(550, HmfLogger.DEBUG_2.intLevel());
    }

    @Test
    public void parseLevelResolvesStandardLevels()
    {
        assertEquals(Level.INFO, HmfLogger.parseLevel("INFO"));
        assertEquals(Level.TRACE, HmfLogger.parseLevel("TRACE"));
    }

    @Test
    public void debug2SitsBetweenDebugAndTrace()
    {
        // finer than DEBUG (shows everything DEBUG does), but coarser than TRACE (excludes third-party TRACE noise)
        assertTrue(HmfLogger.DEBUG_2.isLessSpecificThan(Level.DEBUG));
        assertTrue(HmfLogger.DEBUG_2.isMoreSpecificThan(Level.TRACE));
    }

    @Test
    public void debug2MethodLogsWithoutError()
    {
        // smoke: the extended-logger method resolves and runs with parameterized placeholders
        HmfLogger logger = HmfLogger.getLogger(HmfLoggerTest.class);
        logger.debug2("debug2 smoke {} {}", "a", 1);
    }
}
