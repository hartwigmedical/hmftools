package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UNKNOWN;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class CodingContextTest
{
    @Test
    public void testPrePostCoding()
    {
        // SNVs and MNVs
        int[] exonStarts = { 100, 200, 300, 400, 500, 600 };

        // codons start on at 10, 13, 16 etc
        Integer codingStart = new Integer(350);
        Integer codingEnd = new Integer(425);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 50, codingStart, codingEnd, false, "");

        // upstream
        int pos = 50;
        VariantData var = new VariantData(CHR_1, pos, "A", "C");

        CodingContext codingContext = CodingContext.determineContext(var, transDataPosStrand);

        assertEquals(UPSTREAM, codingContext.RegionType);
        assertEquals(UNKNOWN, codingContext.CodingType);
        assertEquals(PHASE_NONE, codingContext.UpstreamPhase);
        assertEquals(0, codingContext.ExonRank);
        assertEquals(0, codingContext.CodingBase);
        assertEquals(-50, codingContext.NonCodingBaseDistance);

        // 5' UTR exonic
        pos = 125;
        var = new VariantData(CHR_1, pos, "A", "C");
        codingContext = CodingContext.determineContext(var, transDataPosStrand);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_5P, codingContext.CodingType);
        assertEquals(PHASE_NONE, codingContext.UpstreamPhase);
        assertEquals(1, codingContext.ExonRank);
        assertEquals(0, codingContext.CodingBase);
        assertEquals(-126, codingContext.NonCodingBaseDistance);

        // 5' UTR intronic
        pos = 175;
        var = new VariantData(CHR_1, pos, "A", "C");
        codingContext = CodingContext.determineContext(var, transDataPosStrand);

        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(UTR_5P, codingContext.CodingType);
        assertEquals(PHASE_NONE, codingContext.UpstreamPhase);
        assertEquals(1, codingContext.ExonRank);
        assertEquals(0, codingContext.CodingBase);
        assertEquals(-101, codingContext.NonCodingBaseDistance);

        // 5' UTR exonic in same exon as coding begins
        pos = 325;
        var = new VariantData(CHR_1, pos, "A", "C");
        codingContext = CodingContext.determineContext(var, transDataPosStrand);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_5P, codingContext.CodingType);
        assertEquals(PHASE_NONE, codingContext.UpstreamPhase);
        assertEquals(3, codingContext.ExonRank);
        assertEquals(0, codingContext.CodingBase);
        assertEquals(-25, codingContext.NonCodingBaseDistance);

        // 3'UTR exonic in same exon as coding ends
        pos = 430;
        var = new VariantData(CHR_1, pos, "A", "C");
        codingContext = CodingContext.determineContext(var, transDataPosStrand);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_3P, codingContext.CodingType);
        assertEquals(PHASE_NONE, codingContext.UpstreamPhase);
        assertEquals(4, codingContext.ExonRank);
        assertEquals(27, codingContext.CodingBase);
        assertEquals(5, codingContext.NonCodingBaseDistance);

        // 3'UTR intronic
        pos = 555;
        var = new VariantData(CHR_1, pos, "A", "C");
        codingContext = CodingContext.determineContext(var, transDataPosStrand);

        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(UTR_3P, codingContext.CodingType);
        assertEquals(PHASE_NONE, codingContext.UpstreamPhase);
        assertEquals(5, codingContext.ExonRank);
        assertEquals(27, codingContext.CodingBase);
        assertEquals(76, codingContext.NonCodingBaseDistance);

        // 3'UTR exonic
        pos = 625;
        var = new VariantData(CHR_1, pos, "A", "C");
        codingContext = CodingContext.determineContext(var, transDataPosStrand);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_3P, codingContext.CodingType);
        assertEquals(PHASE_NONE, codingContext.UpstreamPhase);
        assertEquals(6, codingContext.ExonRank);
        assertEquals(27, codingContext.CodingBase);
        assertEquals(101, codingContext.NonCodingBaseDistance);
    }

}
