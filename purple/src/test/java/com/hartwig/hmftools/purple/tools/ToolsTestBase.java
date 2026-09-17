package com.hartwig.hmftools.purple.tools;

import java.nio.file.Paths;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency;

public class ToolsTestBase
{
    final String mPurpleDataDir = Paths.get("src", "test", "resources", "tools").toAbsolutePath().toString();
    final String sample1 = "sample1";
    final String sample3 = "sample3";

    RegionAmpDel rad(HumanChromosome chromosome, int start, int stop, AmpDelRegionFrequency.EventType type)
    {
        return new RegionAmpDel(cbr(chromosome, start, stop), type);
    }

    RegionGeneEvents rge(HumanChromosome chromosome, int start, int end, GermlineStatus status)
    {
        return new RegionGeneEvents(gad("G", chromosome, start, end, status));
    }

    RegionGeneEvents rge(String geneName, HumanChromosome chromosome, int start, int end, GermlineStatus status)
    {
        return new RegionGeneEvents(gad(geneName, chromosome, start, end, status));
    }

    GermlineAmpDel gad(String geneName, HumanChromosome chromosome, int start, int end, GermlineStatus status)
    {
        return new GermlineAmpDel(geneName, "TRANSCRIPT1", chromosome.shortName(), "B1", start, end, 0, 0, 0,
                false, null, status, null, 0.0, 0, "", 0, null);
    }

    ChrBaseRegion cbr(HumanChromosome chromosome, int start, int end)
    {
        return new ChrBaseRegion(chromosome.shortName(), start, end);
    }
}
