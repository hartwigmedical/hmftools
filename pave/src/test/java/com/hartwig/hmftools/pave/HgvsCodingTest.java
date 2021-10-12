package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static junit.framework.TestCase.assertEquals;

import org.junit.Test;

public class HgvsCodingTest
{
    @Test
    public void testPointMutations()
    {
        int pos = 100;
        String ref = "A";
        String alt = "C";
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        CodingContext codingContext = new CodingContext();
        codingContext.Strand = POS_STRAND;
        codingContext.RegionType = UPSTREAM;
        codingContext.NearestExonDistance = -50;

        assertEquals("c.-50A>C", HgvsCoding.generate(var, codingContext));

        codingContext.RegionType = INTRONIC;
        codingContext.CodingType = UTR_5P;
        codingContext.NearestExonDistance = -50;

        assertEquals("c.-50A>C", HgvsCoding.generate(var, codingContext));

        codingContext.RegionType = INTRONIC;
        codingContext.CodingType = UTR_3P;
        codingContext.NearestExonDistance = -50;

        assertEquals("c.*50A>C", HgvsCoding.generate(var, codingContext));

        codingContext.RegionType = INTRONIC;
        codingContext.CodingType = CODING;
        codingContext.CodingBase = 10;
        codingContext.NearestExonDistance = -50;

        assertEquals("c.10-50A>C", HgvsCoding.generate(var, codingContext));

        codingContext.NearestExonDistance = 50;

        assertEquals("c.10+50A>C", HgvsCoding.generate(var, codingContext));

        codingContext.RegionType = EXONIC;
        assertEquals("c.10A>C", HgvsCoding.generate(var, codingContext));
    }
}
