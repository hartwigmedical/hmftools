package com.hartwig.hmftools.common.ratio.txt;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.copynumber.cnv.CNVFileLoader;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.ratio.Ratio;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class RatioFileLoaderTest {
    private static final String BASE_PATH = Resources.getResource("ratio").getPath();
    private static final String SAMPLE = "sample";

    @Test
    public void canLoadNormalFile() throws IOException, HartwigException {
        final List<Ratio> ratios = RatioFileLoader.loadNormalRatios(BASE_PATH, SAMPLE);
        assertEquals(1, ratios.size());
    }

    @Test
    public void canLoadTumorFile() throws IOException, HartwigException {
        final List<Ratio> ratios = RatioFileLoader.loadTumorRatios(BASE_PATH, SAMPLE);
        assertEquals(1, ratios.size());
    }
}