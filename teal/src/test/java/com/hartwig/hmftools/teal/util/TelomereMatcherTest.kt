package com.hartwig.hmftools.teal.util

import com.hartwig.hmftools.teal.util.TelomereMatcher.findGTelomereSegment
import org.junit.Before
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.core.config.Configurator
import org.junit.Test

const val EPS = 1e-5

internal class TelomereMatcherTest
{
    @Before
    fun setUp()
    {
        Configurator.setRootLevel(Level.TRACE)
    }

    @Test
    fun testCalcTelomereMatch()
    {
        /*TestCase.assertEquals(1.0, calcGTelomereMatch("TTAGGG"), EPS)
        TestCase.assertEquals(1.0, calcGTelomereMatch("GGTTAGGGT"), EPS)
        TestCase.assertEquals(0.77778, calcGTelomereMatch("GGTAGGT"), EPS)
        TestCase.assertEquals(0.83333, calcGTelomereMatch("TCAGGGGTTAGG"), EPS)*/

        findGTelomereSegment("GTTATTCAACTCCGTGATTTCATTCTATTAGGGTTAGGGTTAGGGTTAGGGTT", 0.9)
    }
}
