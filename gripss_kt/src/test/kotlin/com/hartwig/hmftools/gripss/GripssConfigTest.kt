package com.hartwig.hmftools.gripss

import junit.framework.Assert.assertTrue
import org.junit.Test

class GripssConfigTest {

    @Test
    fun testFilterConfigParametersAllowArguments() {
        val filterOptions = GripssFilterConfig.createOptions()
        for (option in filterOptions.options) {
            assertTrue(option.hasArg())
        }
    }

}