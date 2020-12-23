package com.hartwig.hmftools.common.gene;

import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.calcCodingBases;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.GeneTestUtils;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.CodingBaseData;

import org.junit.Test;

public class GeneTranscriptTest
{
    private static final int TRANS_ID_1 = 1;
    private static final String TRANS_NAME_1 = "TRANS001";
    private static final String GENE_ID_1 = "GENE001";
    @Test
    public void testCodingData()
    {
        Integer codingStart = null;
        Integer codingEnd = null;

        TranscriptData transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, POS_STRAND,
                10, 40, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 1, -1, -1));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 2, -1, -1));

        CodingBaseData cbData = calcCodingBases(transData, 25);

        assertEquals(0, cbData.TotalCodingBases);
        assertEquals(0, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        codingStart = 35;
        codingEnd = 75;

        transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, POS_STRAND,
                10, 80, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 1, PHASE_NONE, PHASE_NONE));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 2, PHASE_NONE, PHASE_0));
        transData.exons().add(new ExonData(TRANS_ID_1, 50, 60, 2, PHASE_0, PHASE_2));
        transData.exons().add(new ExonData(TRANS_ID_1, 70, 80, 2, PHASE_2, PHASE_NONE));

        cbData = calcCodingBases(transData, 25);

        assertEquals(23, cbData.TotalCodingBases);
        assertEquals(0, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        cbData = calcCodingBases(transData, 76);

        assertEquals(cbData.TotalCodingBases, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        cbData = calcCodingBases(transData, 35);

        assertEquals(1, cbData.CodingBases);
        assertEquals(PHASE_1, cbData.Phase);

        cbData = calcCodingBases(transData, 75);

        assertEquals(cbData.TotalCodingBases, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        cbData = calcCodingBases(transData, 45);

        assertEquals(6, cbData.CodingBases);
        assertEquals(PHASE_0, cbData.Phase);

        cbData = calcCodingBases(transData, 65);

        assertEquals(17, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        // test negative strand
        transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, NEG_STRAND,
                10, 80, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 1, PHASE_NONE, PHASE_NONE));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 2, PHASE_2, PHASE_NONE));
        transData.exons().add(new ExonData(TRANS_ID_1, 50, 60, 2, PHASE_0, PHASE_2));
        transData.exons().add(new ExonData(TRANS_ID_1, 70, 80, 2, PHASE_NONE, PHASE_0));

        cbData = calcCodingBases(transData, 25);

        assertEquals(23, cbData.TotalCodingBases);
        assertEquals(cbData.TotalCodingBases, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        cbData = calcCodingBases(transData, 76);

        assertEquals(0, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        cbData = calcCodingBases(transData, 35);

        assertEquals(cbData.TotalCodingBases, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        cbData = calcCodingBases(transData, 75);

        assertEquals(1, cbData.CodingBases);
        assertEquals(PHASE_1, cbData.Phase);

        cbData = calcCodingBases(transData, 45);

        assertEquals(17, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        cbData = calcCodingBases(transData, 65);

        assertEquals(6, cbData.CodingBases);
        assertEquals(PHASE_0, cbData.Phase);

        // test coding start on the first base of the first exon with specific phasing
        codingStart = 10;
        codingEnd = 50;
        transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, POS_STRAND,
                10, 80, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 1, PHASE_2, PHASE_1));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 2, PHASE_1, PHASE_0));
        transData.exons().add(new ExonData(TRANS_ID_1, 50, 60, 2, PHASE_0, PHASE_NONE));

        cbData = calcCodingBases(transData, 10);

        assertEquals(1, cbData.CodingBases);
        assertEquals(PHASE_0, cbData.Phase);

        cbData = calcCodingBases(transData, 20);

        assertEquals(11, cbData.CodingBases);
        assertEquals(PHASE_1, cbData.Phase);

        // and on the negative strand
        codingStart = 30;
        codingEnd = 60;
        transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, NEG_STRAND,
                10, 80, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 1, PHASE_NONE, PHASE_NONE));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 2, PHASE_1, PHASE_0));
        transData.exons().add(new ExonData(TRANS_ID_1, 50, 60, 2, PHASE_2, PHASE_1));

        cbData = calcCodingBases(transData, 60);

        assertEquals(1, cbData.CodingBases);
        assertEquals(PHASE_0, cbData.Phase);

        cbData = calcCodingBases(transData, 40);

        assertEquals(12, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);
    }

    @Test
    public void testTranscriptCreation()
    {
        int[] exonStarts = {10, 30, 50};
        Integer codingStart = 15;
        Integer codingEnd = 55;

        TranscriptData transData = GeneTestUtils.createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, BIOTYPE_PROTEIN_CODING);

        assertEquals(3, transData.exons().size());
        assertEquals(10, transData.exons().get(0).Start);
        assertEquals(20, transData.exons().get(0).End);
        assertEquals(3, transData.exons().get(2).Rank);

        assertEquals(PHASE_NONE, transData.exons().get(0).PhaseStart);
        assertEquals(PHASE_0, transData.exons().get(0).PhaseEnd);
        assertEquals(PHASE_0, transData.exons().get(1).PhaseStart);
        assertEquals(PHASE_2, transData.exons().get(1).PhaseEnd);
        assertEquals(PHASE_2, transData.exons().get(2).PhaseStart);
        assertEquals(PHASE_NONE, transData.exons().get(2).PhaseEnd);

        transData = GeneTestUtils.createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, BIOTYPE_PROTEIN_CODING);

        assertEquals(3, transData.exons().size());
        assertEquals(50, transData.exons().get(2).Start);
        assertEquals(60, transData.exons().get(2).End);
        assertEquals(3, transData.exons().get(0).Rank);

        assertEquals(PHASE_NONE, transData.exons().get(2).PhaseStart);
        assertEquals(PHASE_0, transData.exons().get(2).PhaseEnd);
        assertEquals(PHASE_0, transData.exons().get(1).PhaseStart);
        assertEquals(PHASE_2, transData.exons().get(1).PhaseEnd);
        assertEquals(PHASE_2, transData.exons().get(0).PhaseStart);
        assertEquals(PHASE_NONE, transData.exons().get(0).PhaseEnd);
    }
}
