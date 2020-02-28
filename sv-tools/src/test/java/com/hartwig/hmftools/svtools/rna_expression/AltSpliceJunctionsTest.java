package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.INTRONIC;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NOVEL_3_PRIME;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NOVEL_5_PRIME;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NOVEL_EXON;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.NOVEL_INTRON;
import static com.hartwig.hmftools.svtools.rna_expression.AltSpliceJunctionType.SKIPPED_EXONS;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpressionTest.REF_BASE_STR_1;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpressionTest.createCigar;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpressionTest.createReadRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.junit.Test;

public class AltSpliceJunctionsTest
{
    @Test
    public void testAltSpliceJunctionTypes()
    {
        RnaExpConfig config = new RnaExpConfig();
        config.ExpRateFragmentLengths.add(new int[] { 30, 1 });
        config.ReadLength = 10;

        String chromosome = "1";
        String geneId = "GENE01";

        EnsemblGeneData geneData = new EnsemblGeneData(geneId, geneId, chromosome, (byte) 1, 100, 1500, "");

        int transId1 = 1;
        String transName1 = "TRANS01";

        TranscriptData transData1 = new TranscriptData(transId1, transName1, geneId, true, (byte) 1,
                100, 300, null, null, "");

        transData1.exons().add(new ExonData(transId1, 100, 300, 1, -1, -1));
        transData1.exons().add(new ExonData(transId1, 400, 500, 2, -1, -1));
        transData1.exons().add(new ExonData(transId1, 800, 1000, 3, -1, -1));
        transData1.exons().add(new ExonData(transId1, 1400, 1500, 4, -1, -1));

        int transId2 = 2;
        String transName2 = "TRANS01";

        TranscriptData transData2 = new TranscriptData(transId2, transName2, geneId, true, (byte) 1,
                100, 300, null, null, "");

        transData2.exons().add(new ExonData(transId2, 250, 350, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 400, 500, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 600, 700, 3, -1, -1));

        GeneReadData gene = new GeneReadData(geneData);

        List<TranscriptData> transcripts = Lists.newArrayList(transData1, transData2);
        List<Integer> transIds = Lists.newArrayList(transId1, transId2);

        gene.setTranscripts(transcripts);
        gene.generateExonicRegions();

        AltSpliceJunctionFinder asjFinder = new AltSpliceJunctionFinder(new RnaExpConfig(), null, null);
        asjFinder.setGeneData(gene);

        long[] spliceJunction = new long[SE_PAIR];

        spliceJunction[SE_START] = 300;
        spliceJunction[SE_END] = 400;
        assertTrue(asjFinder.junctionMatchesGene(spliceJunction, transIds));

        spliceJunction[SE_START] = 300;
        spliceJunction[SE_END] = 401;
        assertFalse(asjFinder.junctionMatchesGene(spliceJunction, transIds));

        spliceJunction[SE_START] = 350;
        spliceJunction[SE_END] = 800;
        assertFalse(asjFinder.junctionMatchesGene(spliceJunction, transIds));

        // test out various types of novel junctions:

        // known 5' to novel intronic 3'
        ReadRecord read = createReadRecord(1, chromosome, 491, 559, REF_BASE_STR_1, createCigar(0, 10, 49, 10, 0));

        List<RegionReadData> overlappingRegions = gene.findOverlappingRegions(read);

        read.processOverlappingRegions(overlappingRegions);
        AltSpliceJunction altSJ = asjFinder.createFromRead(read, transIds);

        assertEquals(NOVEL_3_PRIME, altSJ.type());
        assertEquals(0, altSJ.calcNearestExonBoundary(SE_START));
        assertEquals(50, altSJ.calcNearestExonBoundary(SE_END));

        // novel intronic 5' to known 3' prime
        read = createReadRecord(1, chromosome, 361, 409, REF_BASE_STR_1, createCigar(0, 10, 29, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);
        read.processOverlappingRegions(overlappingRegions);
        altSJ = asjFinder.createFromRead(read, transIds);

        assertEquals(NOVEL_5_PRIME, altSJ.type());

        // intronic in both
        read = createReadRecord(1, chromosome, 721, 779, REF_BASE_STR_1, createCigar(0, 10, 39, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);
        assertTrue(overlappingRegions.isEmpty());
        altSJ = asjFinder.createFromRead(read, transIds);

        assertEquals(INTRONIC, altSJ.type());
        assertEquals(30, altSJ.calcNearestExonBoundary(SE_START));
        assertEquals(30, altSJ.calcNearestExonBoundary(SE_END));

        // skipped exons
        read = createReadRecord(1, chromosome, 291, 809, REF_BASE_STR_1, createCigar(0, 10, 499, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);
        read.processOverlappingRegions(overlappingRegions);
        altSJ = asjFinder.createFromRead(read, transIds);

        assertEquals(SKIPPED_EXONS, altSJ.type());
        assertEquals(0, altSJ.calcNearestExonBoundary(SE_START));
        assertEquals(0, altSJ.calcNearestExonBoundary(SE_END));

        // novel intron
        read = createReadRecord(1, chromosome, 841, 959, REF_BASE_STR_1, createCigar(0, 10, 99, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);
        read.processOverlappingRegions(overlappingRegions);
        altSJ = asjFinder.createFromRead(read, transIds);

        assertEquals(NOVEL_INTRON, altSJ.type());
        assertEquals(-150, altSJ.calcNearestExonBoundary(SE_START));
        assertEquals(-150, altSJ.calcNearestExonBoundary(SE_END));

        // fragment reads making a novel exon
        ReadRecord read1 = createReadRecord(1, chromosome, 991, 1209, REF_BASE_STR_1, createCigar(0, 10, 199, 10, 0));
        ReadRecord read2 = createReadRecord(1, chromosome, 1291, 1409, REF_BASE_STR_1, createCigar(0, 10, 99, 10, 0));

        read1.processOverlappingRegions(gene.findOverlappingRegions(read1));
        read2.processOverlappingRegions(gene.findOverlappingRegions(read2));

        transIds = Lists.newArrayList(transId1);
        AltSpliceJunction firstAltSJ = asjFinder.createFromRead(read1, transIds);
        AltSpliceJunction secondAltSJ = asjFinder.createFromRead(read2, transIds);

        asjFinder.checkNovelExon(firstAltSJ, secondAltSJ);

        assertEquals(NOVEL_EXON, firstAltSJ.type());
        assertEquals(NOVEL_EXON, secondAltSJ.type());

    }
}
