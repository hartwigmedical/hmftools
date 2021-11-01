package com.hartwig.hmftools.telo.util;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import org.apache.logging.log4j.Level;
import org.junit.Before;
import org.junit.Test;

public class TelomereMatcherTest
{
    private static final double EPS = 1e-5;

    @Before
    public void setUp()
    {
        org.apache.logging.log4j.core.config.Configurator.setRootLevel(Level.TRACE);
    }

    @Test
    public void testCalcTelomereMatch()
    {
        assertEquals(1.0, TelomereMatcher.calcGTelomereMatch("TTAGGG"), EPS);
        assertEquals(1.0, TelomereMatcher.calcGTelomereMatch("GGTTAGGGT"), EPS);
        assertEquals(0.77778, TelomereMatcher.calcGTelomereMatch("GGTAGGT"), EPS);
        assertEquals(0.83333, TelomereMatcher.calcGTelomereMatch("TCAGGGGTTAGG"), EPS);
    }
}

