package com.hartwig.hmftools.purple;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;

import org.junit.Before;
import org.junit.Test;

public class TargetRegionsCopyNumbersTest extends TargetRegionsTestBase
{
    private final String cobaltDataFilePath = Resources.getResource("cobalt/sample_1.ratios.tsv").getPath();
    private final String purpleDataFilePath = Resources.getResource("purple/segments_1.tsv").getPath();

    @Before
    public void setup()
    {
        super.setup();
    }

    @Test
    public void copyNumbersDataTest() throws IOException
    {
        Map<Chromosome, List<CobaltRatio>> cobaltData = CobaltRatioFile.readWithGender(cobaltDataFilePath, Gender.FEMALE, true);
        List<PurpleCopyNumber> purpleCopyNumbers = PurpleCopyNumberFile.read(purpleDataFilePath);
        TargetRegionsCopyNumbers trc = new TargetRegionsCopyNumbers(targetRegionsData, cobaltData, purpleCopyNumbers, refGenomeVersion);

        List<TargetRegionsCopyNumber> copyNumberData = trc.copyNumbersData();
        assertEquals(4, copyNumberData.size());

        TargetRegionsCopyNumber cn0 = copyNumberData.get(0);
        assertEquals(26694001,cn0.mCobaltRatio().position());
        assertEquals(3,cn0.mOverlappingRegions().size());
        assertEquals(26694021,cn0.mOverlappingRegions().get(0).start());
        assertEquals("ARID1A_0",cn0.mOverlappingRegions().get(0).mTag);
        assertEquals("ARID1A_1",cn0.mOverlappingRegions().get(1).mTag);
        assertEquals("ARID1A_2",cn0.mOverlappingRegions().get(2).mTag);
        assertEquals(1, cn0.mPurpleCopyNumber().start());
        assertEquals(63461950, cn0.mPurpleCopyNumber().end());

        TargetRegionsCopyNumber cn2 = copyNumberData.get(2);
        assertEquals(55020001,cn2.mCobaltRatio().position());
        assertEquals(1,cn2.mOverlappingRegions().size());
        assertEquals(55020021,cn2.mOverlappingRegions().get(0).start());
        assertEquals("EGFR_1",cn2.mOverlappingRegions().get(0).mTag);
        assertEquals(38402474, cn2.mPurpleCopyNumber().start());
        assertEquals(59554330, cn2.mPurpleCopyNumber().end());

        TargetRegionsCopyNumber cn3 = copyNumberData.get(3);
        assertEquals(55021001,cn3.mCobaltRatio().position());
        assertEquals(1,cn3.mOverlappingRegions().size());
        assertEquals(55020021,cn3.mOverlappingRegions().get(0).start());
        assertEquals(55021520,cn3.mOverlappingRegions().get(0).end());
        assertEquals("EGFR_1",cn3.mOverlappingRegions().get(0).mTag);
    }
}
