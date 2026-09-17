package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._20;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._21;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._22;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.file.Paths;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;

import org.junit.Test;

public class RegionEventsExtractorTest
{
    private final String mPurpleDataDir = Paths.get("src", "test", "resources", "tools").toAbsolutePath().toString();

    @Test
    public void runFor1Sample() throws IOException
    {
        final String sample3 = "sample3";
        RegionEventExtractor extractor = new RegionEventExtractor(dataFile(sample3));
        ListMultimap<HumanChromosome, RegionGeneEvents> extracted = extractor.events();
        assertEquals(3, extracted.asMap().size());
        check(extracted.get(_20).get(0), _20, 58692001, 58693000);
        check(extracted.get(_20).get(1), _20, 58981001, 58982000);
        check(extracted.get(_21).get(0), _21, 7929001, 8550172);
        check(extracted.get(_22).get(0), _22, 10699413, 12694000);
    }

    private void check(RegionGeneEvents rge, HumanChromosome chromosome, int start, int end)
    {
        assertEquals(chromosome, rge.chromosome());
        assertEquals(start, rge.region().start());
        assertEquals(end, rge.region().end());
    }

    private String dataFile(String sample)
    {
        return GermlineAmpDel.generateFilename(mPurpleDataDir, sample);
    }
}
