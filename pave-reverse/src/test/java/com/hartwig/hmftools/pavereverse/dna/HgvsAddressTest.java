package com.hartwig.hmftools.pavereverse.dna;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;
import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

import org.junit.Test;

public class HgvsAddressTest extends ReversePaveTestBase
{
    int[] exonStarts = { 10, 30, 50 };
    int codingStart = 15;
    int codingEnd = 55;
    GeneData geneData = createEnsemblGeneData("id132", "BLAH", "1",  POS_STRAND, 5, 105);
    TranscriptData transcript = createTransExons(geneData.GeneId, 123, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "whatever");
    GeneTranscript gt = new GeneTranscript(geneData, transcript);
    GeneData geneDataRS = createEnsemblGeneData("id132", "BLAH", "1",  NEG_STRAND, 5, 105);
    TranscriptData transcriptRS = createTransExons(geneDataRS.GeneId, 123, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "whatever");
    GeneTranscript gtRS = new GeneTranscript(geneDataRS, transcriptRS);
    // 12345678901234567890123456789012345678901234567890123456789012345678901
    // 1        10        20        30        40        50        60        70
    // ____-----+++++******---------***********---------******+++++-----_____

    int[] exonStarts2 = { 10, 30, 50, 70, 90 };
    int codingStart2 = 35;
    int codingEnd2 = 74;
    GeneData geneData2 = createEnsemblGeneData("id132", "BLAH", "1",  POS_STRAND, 5, 105);
    TranscriptData transcript2 = createTransExons(geneData2.GeneId, 123, POS_STRAND, exonStarts2, 9, codingStart2, codingEnd2, false, "whatever");
    GeneTranscript gt2 = new GeneTranscript(geneData2, transcript2);
    GeneData geneData2RS = createEnsemblGeneData("id132", "BLAH", "1",  NEG_STRAND, 5, 105);
    TranscriptData transcript2RS = createTransExons(geneData2.GeneId, 123, NEG_STRAND, exonStarts2, 9, codingStart2, codingEnd2, false, "whatever");
    GeneTranscript gt2RS = new GeneTranscript(geneData2RS, transcript2RS);
    // 12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    // 1        10        20        30        40        50        60        70        80        90        100       110
    // ____-----++++++++++----------+++++*****----------**********----------*****+++++----------++++++++++------_____

    int[] exonStarts3 = { 10, 30, 50, 70, 90 };
    int codingStart3 = 55;
    int codingEnd3 = 94;
    GeneData geneData3 = createEnsemblGeneData("id132", "BLAH", "1",  POS_STRAND, 5, 105);
    TranscriptData transcript3 = createTransExons(geneData3.GeneId, 123, POS_STRAND, exonStarts3, 9, codingStart3, codingEnd3, false, "whatever");
    GeneTranscript gt3 = new GeneTranscript(geneData3, transcript3);
    GeneData geneData3RS = createEnsemblGeneData("id132", "BLAH", "1",  NEG_STRAND, 5, 105);
    TranscriptData transcript3RS = createTransExons(geneData3RS.GeneId, 123, NEG_STRAND, exonStarts3, 9, codingStart3, codingEnd3, false, "whatever");
    GeneTranscript gt3RS = new GeneTranscript(geneData3RS, transcript3RS);
    // 12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    // 1        10        20        30        40        50        60        70        80        90        100       110
    // ____-----++++++++++----------++++++++++----------+++++*****----------**********----------*****+++++------_____

    @Test
    public void inIntronBeforeExonReverseStrand()
    {
        assertEquals(41, new InIntronBeforeExon(7, 1).toStrandLocation(gtRS));
        assertEquals(42, new InIntronBeforeExon(7, 2).toStrandLocation(gtRS));
        assertEquals(21, new InIntronBeforeExon(18, 1).toStrandLocation(gtRS));
    }

    @Test
    public void inIntronBeforeExon()
    {
        assertEquals(30 - 1, new InIntronBeforeExon(7, 1).toStrandLocation(gt));
        assertEquals(30 - 2, new InIntronBeforeExon(7, 2).toStrandLocation(gt));
        assertEquals(50 - 1, new InIntronBeforeExon(18, 1).toStrandLocation(gt));
        assertEquals(50 - 3, new InIntronBeforeExon(18, 3).toStrandLocation(gt));
    }

    @Test
    public void inIntronUpstreamOfCodingStartReverseStrand()
    {
        assertEquals(80, new InIntronUpstreamOfCodingStart(-5, -1).toStrandLocation(gt2RS));
        assertEquals(81, new InIntronUpstreamOfCodingStart(-5, -2).toStrandLocation(gt2RS));
        assertEquals(100, new InIntronUpstreamOfCodingStart(-15, -1).toStrandLocation(gt2RS));
    }

    @Test
    public void inIntronUpstreamOfCodingStart()
    {
        assertEquals(49, new InIntronUpstreamOfCodingStart(-5, -1).toStrandLocation(gt3));
        assertEquals(48, new InIntronUpstreamOfCodingStart(-5, -2).toStrandLocation(gt3));
        assertEquals(40, new InIntronUpstreamOfCodingStart(-6, 1).toStrandLocation(gt3));
        assertEquals(41, new InIntronUpstreamOfCodingStart(-6, 2).toStrandLocation(gt3));
        assertEquals(29, new InIntronUpstreamOfCodingStart(-15, -1).toStrandLocation(gt3));
        assertEquals(20, new InIntronUpstreamOfCodingStart(-16, 1).toStrandLocation(gt3));
        assertEquals(9, new InIntronUpstreamOfCodingStart(-25, -1).toStrandLocation(gt3));
    }

    @Test
    public void inIntronDownstreamOfCodingEnd()
    {
        assertEquals(80, new InIntronDownstreamOfCodingEnd(5, 1).toStrandLocation(gt2));
        assertEquals(81, new InIntronDownstreamOfCodingEnd(5, 2).toStrandLocation(gt2));
        assertEquals(89, new InIntronDownstreamOfCodingEnd(6, -1).toStrandLocation(gt2));
        assertEquals(88, new InIntronDownstreamOfCodingEnd(6, -2).toStrandLocation(gt2));
        assertEquals(100, new InIntronDownstreamOfCodingEnd(15, 1).toStrandLocation(gt2));
    }

    @Test
    public void inIntronDownstreamOfCodingEndReverseStrand()
    {
        assertEquals(29, new InIntronDownstreamOfCodingEnd(5, 1).toStrandLocation(gt2RS));
        assertEquals(28, new InIntronDownstreamOfCodingEnd(5, 2).toStrandLocation(gt2RS));
        assertEquals(9, new InIntronDownstreamOfCodingEnd(15, 1).toStrandLocation(gt2RS));
    }

    @Test
    public void inIntronAfterExon()
    {
        assertEquals(21, new InIntronAfterExon(6, 1).toStrandLocation(gt));
        assertEquals(22, new InIntronAfterExon(6, 2).toStrandLocation(gt));
        assertEquals(30 + 10 + 1, new InIntronAfterExon(17, 1).toStrandLocation(gt));
        assertEquals(30 + 10 + 3, new InIntronAfterExon(17, 3).toStrandLocation(gt));
    }

    @Test
    public void inIntronAfterExonReverseStrand()
    {
        assertEquals(49, new InIntronAfterExon(6, 1).toStrandLocation(gtRS));
        assertEquals(48, new InIntronAfterExon(6, 2).toStrandLocation(gtRS));
        assertEquals(29, new InIntronAfterExon(17, 1).toStrandLocation(gtRS));
        assertEquals(27, new InIntronAfterExon(17, 3).toStrandLocation(gtRS));
    }

    @Test
    public void inExon()
    {
        assertEquals(15, new InExon(1).toStrandLocation(gt));
        assertEquals(20, new InExon(6).toStrandLocation(gt));
        assertEquals(30, new InExon(7).toStrandLocation(gt));
        assertEquals(30 + 10, new InExon(17).toStrandLocation(gt));
        assertEquals(50, new InExon(18).toStrandLocation(gt));
        assertEquals(55, new InExon(23).toStrandLocation(gt));
    }

    @Test
    public void inExonReverseStrand()
    {
        assertEquals(55, new InExon(1).toStrandLocation(gtRS));
        assertEquals(40, new InExon(7).toStrandLocation(gtRS));
    }

    @Test
    public void inExonDownstreamOfCodingEnd()
    {
        assertEquals(56, new InExonDownstreamOfCodingEnd(1).toStrandLocation(gt));
        assertEquals(57, new InExonDownstreamOfCodingEnd(2).toStrandLocation(gt));
    }

    @Test
    public void inExonDownstreamOfCodingEndReverseStrand()
    {
        assertEquals(14, new InExonDownstreamOfCodingEnd(1).toStrandLocation(gtRS));
        assertEquals(13, new InExonDownstreamOfCodingEnd(2).toStrandLocation(gtRS));
    }

    @Test
    public void inExonUpstreamOfCodingStart()
    {
        assertEquals(14, new InExonUpstreamOfCodingStart(-1).toStrandLocation(gt));
        assertEquals(13, new InExonUpstreamOfCodingStart(-2).toStrandLocation(gt));
        assertEquals(10, new InExonUpstreamOfCodingStart(-5).toStrandLocation(gt));

        // This is before the exon, but we need this for some variants (eg TERT)
        assertEquals(9, new InExonUpstreamOfCodingStart(-6).toStrandLocation(gt));
        assertEquals(8, new InExonUpstreamOfCodingStart(-7).toStrandLocation(gt));
    }

    @Test
    public void inExonUpstreamOfCodingStartReverseStrand()
    {
        assertEquals(56, new InExonUpstreamOfCodingStart(-1).toStrandLocation(gtRS));
        assertEquals(57, new InExonUpstreamOfCodingStart(-2).toStrandLocation(gtRS));
        assertEquals(60, new InExonUpstreamOfCodingStart(-5).toStrandLocation(gtRS));

        // This is before the exon, but we need this for some variants (eg TERT)
        assertEquals(61, new InExonUpstreamOfCodingStart(-6).toStrandLocation(gtRS));
        assertEquals(62, new InExonUpstreamOfCodingStart(-7).toStrandLocation(gtRS));
    }

}
