package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.createBreakendTranscriptData;
import static com.hartwig.hmftools.common.gene.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.gene.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.gene.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.gene.GeneTestUtils.createGeneAnnotation;
import static com.hartwig.hmftools.common.gene.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.fusion.FusionCommon.DEFAULT_PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.linx.fusion.rna.RnaFusionAnnotator;

import org.junit.Test;

public class RnaFusionTest
{
    @Test
    public void testRnaMatching()
    {
        EnsemblDataCache geneTransCache = new EnsemblDataCache("", RefGenomeVersion.V37);

        String geneName = "GENE1";
        String geneId = "ENSG0001";
        String chromosome = "1";

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, 1, 10000, 20000));

        // one on the negative strand
        String geneName2 = "GENE2";
        String geneId2 = "ENSG0003";
        geneList.add(createEnsemblGeneData(geneId2, geneName2, chromosome, -1, 10000, 20000));

        addGeneData(geneTransCache, chromosome, geneList);

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;
        byte strand = 1;

        int[] exonStarts = new int[]{10500, 11500, 12500, 13500};
        int codingStart = 10500;
        int codingEnd = 13599;

        TranscriptData transData = createTransExons(geneId, transId++, strand, exonStarts, 100, codingStart, codingEnd, true, "");
        String transName = transData.TransName;
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId, transDataList);

        transDataList = Lists.newArrayList();

        strand = -1;

        transData = createTransExons(geneId, transId++, strand, exonStarts, 100, codingStart, codingEnd, true, "");
        String transName2 = transData.TransName;
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId2, transDataList);

        PRE_GENE_PROMOTOR_DISTANCE = DEFAULT_PRE_GENE_PROMOTOR_DISTANCE;
        RnaFusionAnnotator rnaAnnotator = new RnaFusionAnnotator(geneTransCache);

        // test positive strand
        int svPos1 = 12700;
        BreakendGeneData geneAnnot1 = createGeneAnnotation(0, true, geneName, geneId, POS_STRAND, chromosome, svPos1, 1);

        transData = geneTransCache.getTranscriptData(geneId, transName);
        assertEquals(4, transData.exons().size());

        BreakendTransData trans = createBreakendTranscriptData(transData, geneAnnot1.position(), geneAnnot1);

        assertNotNull(trans);

        // test upstream scenarios
        int rnaPosition = 12600;
        boolean isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans, true, geneAnnot1.position(), rnaPosition, true);

        assertTrue(isValid);

        // after the next splice site
        svPos1 = 13500;
        geneAnnot1.setPositionalData(chromosome, svPos1, (byte)1);

        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans, true, geneAnnot1.position(), rnaPosition, true);

        assertFalse(isValid);

        // test non-exact RNA boundary
        rnaPosition = 12550;
        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans, true, geneAnnot1.position(), rnaPosition, false);

        assertFalse(isValid);

        rnaPosition = 12700;
        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans, true, geneAnnot1.position(), rnaPosition, false);

        assertFalse(isValid);

        // test downstream

        // exact base at 2nd exon
        svPos1 = 100; // pre promotor
        geneAnnot1.setPositionalData(chromosome, svPos1, (byte)1);

        rnaPosition = 11500;
        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans, false, geneAnnot1.position(), rnaPosition, true);

        assertTrue(isValid);

        // before previous splice acceptor
        svPos1 = 12400;
        geneAnnot1.setPositionalData(chromosome, svPos1, (byte)1);

        rnaPosition = 13500;
        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans, false, geneAnnot1.position(), rnaPosition, true);

        assertFalse(isValid);

        svPos1 = 13000;
        geneAnnot1.setPositionalData(chromosome, svPos1, (byte)1);

        // valid position
        rnaPosition = 13500;
        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans, false, geneAnnot1.position(), rnaPosition, true);

        assertTrue(isValid);


        // now test the negative strand

        int svPos2 = 12700;
        BreakendGeneData geneAnnot2 = createGeneAnnotation(1, true, geneName2, geneId2, NEG_STRAND, chromosome, svPos2, -1);

        transData = geneTransCache.getTranscriptData(geneId2, transName2);
        assertEquals(4, transData.exons().size());

        BreakendTransData trans2 = createBreakendTranscriptData(transData, geneAnnot2.position(), geneAnnot2);

        assertNotNull(trans2);

        // upstream

        rnaPosition = 11500; // 3rd exon end

        svPos2 = 11600;
        geneAnnot2.setPositionalData(chromosome, svPos2, (byte)1);

        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans2, true, geneAnnot2.position(), rnaPosition, true);

        assertTrue(isValid);

        // test downstream

        rnaPosition = 11600; // 3rd exon start

        svPos2 = 11700;
        geneAnnot2.setPositionalData(chromosome, svPos2, (byte)1);

        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans2, false, geneAnnot2.position(), rnaPosition, true);

        assertTrue(isValid);

        // before prev splice acceptor is invali
        svPos2 = 12700;
        geneAnnot2.setPositionalData(chromosome, svPos2, (byte)1);

        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans2, false, geneAnnot2.position(), rnaPosition, true);

        assertFalse(isValid);

        // invalid too far upstream of promotor
        rnaPosition = 12600; // 3rd exon start

        svPos2 = 130000;
        geneAnnot2.setPositionalData(chromosome, svPos2, (byte)1);

        isValid = rnaAnnotator.isTranscriptBreakendViableForRnaBoundary(
                trans2, false, geneAnnot2.position(), rnaPosition, true);

        assertFalse(isValid);

    }
}
