package com.hartwig.hmftools.svanalysis.fusion;

import static com.hartwig.hmftools.svanalysis.analyser.com.hartwig.hmftools.svanalysis.gene.GeneTestUtils.createGeneAnnotation;
import static com.hartwig.hmftools.svanalysis.fusion.FusionFinder.checkFusionLogic;

import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.junit.Test;

public class FusionRulesTest
{
    @Test
    public void testFusionTypes()
    {
        String geneName = "GENE1";
        String geneId = "ENSG0001";
        String chromosome = "1";

        // SV breakend positions won't impact fusion determination since transcripts are created manually
        GeneAnnotation gene1 = createGeneAnnotation(0, true, geneName, geneId, 1, chromosome, 0, 1);

        // one on the negative strand
        String geneName2 = "GENE2";
        String geneId2 = "ENSG0003";
        String chromosome2 = "1";

        GeneAnnotation gene2 = createGeneAnnotation(0, false, geneName2, geneId2, -1, chromosome2, 0, 1);

        String transName1 = "ENST0001";
        int transId1 = 1;

        /* invalid fusions:
            - up or down post-coding
            - down non-coding
            - up promoter
            - down single exon
            - up not disruptive
            - up 5'UTR to down coding
            - up coding to down not-coding
            - up coding exonic to not down coding exonic
            - up non-coding to down coding
            - unmatched phasing - exact or not
         */


        // non-coding combos
        Long codingStart = null;
        Long codingEnd = null;

        Transcript trans1 = new Transcript(gene1, transId1, transName1, 2, -1, 3, -1,
                0, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        String transName2 = "ENST0002";
        int transId2 = 2;

        codingStart = new Long(100);
        codingEnd = new Long(200);
        Transcript trans2 = new Transcript(gene2, transId2, transName2, 2, -1, 3, -1,
                10, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        // up non-coding
        assertTrue(trans1.nonCoding());
        assertTrue(trans2.isCoding());
        assertTrue(checkFusionLogic(trans1, trans2) == null);

        codingStart = null;
        codingEnd = null;
        trans2 = new Transcript(gene2, transId2, transName2, 2, -1, 3, -1,
                0, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        codingStart = new Long(100);
        codingEnd = new Long(200);
        trans1 = new Transcript(gene1, transId1, transName1, 2, -1, 3, -1,
                10, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        // down non-coding
        assertTrue(trans1.isCoding());
        assertTrue(trans2.nonCoding());
        assertTrue(checkFusionLogic(trans1, trans2) == null);

        // up / down post coding
        trans2 = new Transcript(gene2, transId2, transName2, 2, -1, 3, -1,
                10, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        trans1 = new Transcript(gene1, transId1, transName1, 2, -1, 3, -1,
                getCodingBases(codingStart, codingEnd), getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        assertTrue(trans1.postCoding());
        assertTrue(trans2.isCoding());
        assertTrue(checkFusionLogic(trans1, trans2) == null);

        // up promotor
        trans1 = new Transcript(gene1, transId1, transName1, 0, -1, 1, -1,
                0, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        assertTrue(trans1.isPromoter());
        assertTrue(checkFusionLogic(trans1, trans2) == null);

        // down single exon
        trans1 = new Transcript(gene1, transId1, transName1, 2, -1, 3, -1,
                0, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        trans2 = new Transcript(gene2, transId2, transName2, 2, -1, 3, -1,
                10, getCodingBases(codingStart, codingEnd),1, true, 50, 250, codingStart, codingEnd);

        assertTrue(trans1.preCoding());
        assertTrue(checkFusionLogic(trans1, trans2) == null);

        // up not disruptive
        trans1.setIsDisruptive(false);

        trans2 = new Transcript(gene2, transId2, transName2, 2, -1, 3, -1,
                10, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        assertTrue(checkFusionLogic(trans1, trans2) == null);

        // up 5'UTR to down coding
        trans1.setIsDisruptive(true);

        // up coding exonic to not down coding exonic
        trans1 = new Transcript(gene1, transId1, transName1, 2, -1, 2, -1,
                0, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        trans2 = new Transcript(gene2, transId2, transName2, 2, -1, 3, -1,
                0, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        assertTrue(trans1.isExonic());
        assertTrue(trans2.isIntronic());
        assertTrue(checkFusionLogic(trans1, trans2) == null);

        // unmatched phasing intronic
        trans1 = new Transcript(gene1, transId1, transName1, 2, -1, 3, -1,
                10, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        trans2 = new Transcript(gene2, transId2, transName2, 2, 0, 3, 0,
                10, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        assertTrue(trans1.isCoding());
        assertTrue(trans2.isCoding());
        assertTrue(trans1.isIntronic());
        assertTrue(trans2.isIntronic());
        assertTrue(checkFusionLogic(trans1, trans2) == null);

        // unmatched phasing exon to exon
        trans1 = new Transcript(gene1, transId1, transName1, 2, 2, 2, 1,
                10, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        trans2 = new Transcript(gene2, transId2, transName2, 2, 1, 2, 0,
                11, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        assertTrue(trans1.isExonic());
        assertTrue(trans2.isExonic());
        assertTrue(checkFusionLogic(trans1, trans2) == null);

        // valid exonic fusion
        trans1 = new Transcript(gene1, transId1, transName1, 2, 2, 2, 1,
                11, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        assertTrue(checkFusionLogic(trans1, trans2) != null);

        // valid intronic fusion
        trans1 = new Transcript(gene1, transId1, transName1, 2, 1, 3, 2,
                10, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        trans2 = new Transcript(gene2, transId2, transName2, 2, 0, 3, 1,
                10, getCodingBases(codingStart, codingEnd),4, true, 50, 250, codingStart, codingEnd);

        assertTrue(checkFusionLogic(trans1, trans2) != null);

    }

    private static int getCodingBases(final Long start, final Long end)
    {
        if(start != null && end != null)
            return (int)(end - start);
        return 0;
    }


}
