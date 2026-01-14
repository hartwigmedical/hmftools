package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.AMP;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.DEL;

import java.io.File;
import java.nio.file.Paths;

import org.junit.Assert;
import org.junit.Test;

public class GermlineAmpDelFrequencyCacheTest
{
    private final String GermlineDirectory = Paths.get("src", "test", "resources", "germline").toAbsolutePath().toString();

    @Test
    public void loadOldFormat()
    {
        File sample1 = new File(GermlineDirectory, "sample1.germline.frequencies.csv");
        GermlineAmpDelFrequencyCache cache = new GermlineAmpDelFrequencyCache(sample1.getAbsolutePath());
        Assert.assertEquals(17, cache.getRegionFrequency("chr9", 214001, 216000, 0, DEL));
        Assert.assertEquals(0, cache.getRegionFrequency("chr9", 214001, 216000, 0, AMP));
    }

    @Test
    public void loadNewFormat()
    {
        File sample2 = new File(GermlineDirectory, "sample2.germline.frequencies.csv");
        GermlineAmpDelFrequencyCache cache = new GermlineAmpDelFrequencyCache(sample2.getAbsolutePath());
        Assert.assertEquals(444, cache.getRegionFrequency("chr9", 400001, 401000, 0, DEL));
        Assert.assertEquals(6, cache.getRegionFrequency("chr9", 400001, 401000, 0, AMP));
    }
}
