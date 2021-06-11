package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MNV;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_THREADS;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SageConfigTest
{
    @Test
    public void testBqrFile()
    {
        SageConfig config = new SageConfig();
        assertEquals("SAMPLE.sage.bqr.tsv", config.baseQualityRecalibrationFile("SAMPLE"));

        config = new SageConfig();
        assertEquals("SAMPLE.sage.bqr.tsv", config.baseQualityRecalibrationFile("SAMPLE"));
    }

}
