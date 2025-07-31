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

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

public class TargetRegionsCopyNumbersTest extends TargetRegionsTestBase
{
    private final String cobaltDataFilePath = Resources.getResource("cobalt/sample_1.ratios.tsv").getPath();

    @Before
    public void setup()
    {
        super.setup();
    }

    @Ignore
    @Test
    public void copyNumbersDataTest() throws IOException
    {
        Map<Chromosome, List<CobaltRatio>> cobaltData = CobaltRatioFile.readWithGender(cobaltDataFilePath, Gender.FEMALE, true);
        TargetRegionsCopyNumbers trc = new TargetRegionsCopyNumbers(targetRegionsData, cobaltData, refGenomeVersion);

        List<TargetRegionsCopyNumber> copyNumberData = trc.copyNumbersData();
        assertEquals(3, copyNumberData.size());
    }
}
