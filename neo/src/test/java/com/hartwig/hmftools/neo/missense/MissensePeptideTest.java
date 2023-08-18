package com.hartwig.hmftools.neo.missense;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.GeneTestUtils;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class MissensePeptideTest
{
    @Test
    public void testMissensePepetides()
    {
        MissenseConfig config = new MissenseConfig(3, 0);
        MockRefGenome refGenome = new MockRefGenome();
        MissenseCalcs calcs = new MissenseCalcs(config, refGenome);

        GeneData testGene = new GeneData(GENE_ID_1, GENE_NAME_1, CHR_1, NEG_STRAND, 10, 80, "");

        int codingStart = 15;
        int codingEnd = 75;
        TranscriptData transDataNeg = GeneTestUtils.createTransExons(
                testGene.GeneId, TRANS_ID_1, NEG_STRAND,
                new int[] {10, 30, 50, 70}, 10, codingStart, codingEnd, true, "");

        calcs.processTranscript(testGene, transDataNeg);

    }
}
