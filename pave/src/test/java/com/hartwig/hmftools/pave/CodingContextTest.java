package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UNKNOWN;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.gene.TranscriptData;

import org.junit.Test;

public class CodingContextTest
{
    private static CodingContext determineContext(final VariantData variant, final TranscriptData transData)
    {
        CodingContext codingContext = new CodingContext();
        CodingUtils.determineContext(variant, transData, codingContext);
        return codingContext;
    }

    private VariantData createVariant(int pos)
    {
        return new VariantData(CHR_1, pos, "A", "C");
    }

    @Test
    public void testNonCoding()
    {
        int[] exonStarts = { 10, 30, 50 };

        TranscriptData transDataPos = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, null, null, false, "");

        // intronic
        int pos = 22;
        VariantData var = createVariant(pos);

        CodingContext codingContext = determineContext(var, transDataPos);

        assertEquals(NON_CODING, codingContext.CodingType);
        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(1, codingContext.ExonRank);
        assertEquals(11, codingContext.CodingBase); // exonic bases since start
        assertEquals(2, codingContext.NearestExonDistance); // 2 bases into next intron

        pos = 47;
        var = createVariant(pos);

        codingContext = determineContext(var, transDataPos);

        assertEquals(NON_CODING, codingContext.CodingType);
        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(3, codingContext.ExonRank);
        assertEquals(23, codingContext.CodingBase); // exonic bases since start
        assertEquals(-3, codingContext.NearestExonDistance); // 2 bases into next intron

        // exonic
        pos = 35;
        var = createVariant(pos);

        codingContext = determineContext(var, transDataPos);

        assertEquals(NON_CODING, codingContext.CodingType);
        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(2, codingContext.ExonRank);
        assertEquals(17, codingContext.CodingBase); // exonic bases since start
        assertEquals(0, codingContext.NearestExonDistance); // 2 bases into next intron

        // negative strand
        TranscriptData transDataNeg = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, null, null, false, "");

        // intronic
        pos = 48;
        var = createVariant(pos);

        codingContext = determineContext(var, transDataNeg);

        assertEquals(NON_CODING, codingContext.CodingType);
        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(1, codingContext.ExonRank);
        assertEquals(11, codingContext.CodingBase); // exonic bases since start
        assertEquals(2, codingContext.NearestExonDistance); // 2 bases into next intron

        pos = 23;
        var = createVariant(pos);

        codingContext = determineContext(var, transDataNeg);

        assertEquals(NON_CODING, codingContext.CodingType);
        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(3, codingContext.ExonRank);
        assertEquals(23, codingContext.CodingBase); // exonic bases since start
        assertEquals(-3, codingContext.NearestExonDistance); // 2 bases into next intron

        // exonic
        pos = 35;
        var = createVariant(pos);

        codingContext = determineContext(var, transDataNeg);

        assertEquals(NON_CODING, codingContext.CodingType);
        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(2, codingContext.ExonRank);
        assertEquals(17, codingContext.CodingBase); // exonic bases since start
        assertEquals(0, codingContext.NearestExonDistance); // 2 bases into next intron
    }

    @Test
    public void testPosStrand()
    {
        int[] exonStarts = { 10, 30, 50, 70, 90 };

        Integer codingStart = new Integer(35);
        Integer codingEnd = new Integer(55);

        TranscriptData transDataPos = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        // upstream
        int pos = 5;
        VariantData var = createVariant(pos);

        CodingContext codingContext = determineContext(var, transDataPos);

        assertEquals(UPSTREAM, codingContext.RegionType);
        assertEquals(UNKNOWN, codingContext.CodingType);

        // 5' UTR exonic
        pos = 15;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataPos);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_5P, codingContext.CodingType);
        assertEquals(1, codingContext.ExonRank);
        assertEquals(11, codingContext.CodingBase);
        assertEquals(0, codingContext.NearestExonDistance);

        // 5' UTR intronic
        pos = 22;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataPos);

        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(UTR_5P, codingContext.CodingType);
        assertEquals(1, codingContext.ExonRank);
        assertEquals(6, codingContext.CodingBase);
        assertEquals(2, codingContext.NearestExonDistance);

        // 5' UTR exonic in same exon as coding begins
        pos = 32;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataPos);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_5P, codingContext.CodingType);
        assertEquals(2, codingContext.ExonRank);
        assertEquals(3, codingContext.CodingBase);
        assertEquals(0, codingContext.NearestExonDistance);

        // 3'UTR exonic in same exon as coding ends
        pos = 59;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataPos);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_3P, codingContext.CodingType);
        assertEquals(3, codingContext.ExonRank);
        assertEquals(4, codingContext.CodingBase);
        assertEquals(0, codingContext.NearestExonDistance);

        // 3'UTR intronic
        pos = 82;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataPos);

        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(UTR_3P, codingContext.CodingType);
        assertEquals(4, codingContext.ExonRank);
        assertEquals(16, codingContext.CodingBase);
        assertEquals(2, codingContext.NearestExonDistance);

        // 3'UTR exonic
        pos = 95;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataPos);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_3P, codingContext.CodingType);
        assertEquals(5, codingContext.ExonRank);
        assertEquals(22, codingContext.CodingBase);
        assertEquals(0, codingContext.NearestExonDistance);
    }

    @Test
    public void testNegStrand()
    {
        int[] exonStarts = { 10, 30, 50, 70, 90 };

        Integer codingStart = new Integer(35);
        Integer codingEnd = new Integer(55);

        TranscriptData transDataNeg = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        // upstream
        int pos = 105;
        VariantData var = createVariant(pos);

        CodingContext codingContext = determineContext(var, transDataNeg);

        assertEquals(UPSTREAM, codingContext.RegionType);
        assertEquals(UNKNOWN, codingContext.CodingType);

        // 5' UTR exonic
        pos = 95;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataNeg);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_5P, codingContext.CodingType);
        assertEquals(1, codingContext.ExonRank);
        assertEquals(22, codingContext.CodingBase);
        assertEquals(0, codingContext.NearestExonDistance);

        // 5' UTR intronic
        pos = 88;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataNeg);

        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(UTR_5P, codingContext.CodingType);
        assertEquals(1, codingContext.ExonRank);
        assertEquals(17, codingContext.CodingBase);
        assertEquals(2, codingContext.NearestExonDistance);

        // 5' UTR exonic in same exon as coding begins
        pos = 58;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataNeg);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_5P, codingContext.CodingType);
        assertEquals(3, codingContext.ExonRank);
        assertEquals(3, codingContext.CodingBase);
        assertEquals(0, codingContext.NearestExonDistance);

        // 3'UTR exonic in same exon as coding ends
        pos = 30;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataNeg);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_3P, codingContext.CodingType);
        assertEquals(4, codingContext.ExonRank);
        assertEquals(5, codingContext.CodingBase);
        assertEquals(0, codingContext.NearestExonDistance);

        // 3'UTR intronic
        pos = 22;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataNeg);

        assertEquals(INTRONIC, codingContext.RegionType);
        assertEquals(UTR_3P, codingContext.CodingType);
        assertEquals(5, codingContext.ExonRank);
        assertEquals(6, codingContext.CodingBase);
        assertEquals(-2, codingContext.NearestExonDistance);

        // 3'UTR exonic
        pos = 15;
        var = createVariant(pos);
        codingContext = determineContext(var, transDataNeg);

        assertEquals(EXONIC, codingContext.RegionType);
        assertEquals(UTR_3P, codingContext.CodingType);
        assertEquals(5, codingContext.ExonRank);
        assertEquals(11, codingContext.CodingBase);
        assertEquals(0, codingContext.NearestExonDistance);
    }

}
