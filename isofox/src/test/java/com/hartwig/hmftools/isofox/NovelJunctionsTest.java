package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createIsofoxConfig;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.TestUtils.createSupplementaryReadPair;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionContext.SPLICE_JUNC;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.CIRCULAR;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.INTRONIC;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.MIXED_TRANS;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_3_PRIME;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_5_PRIME;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_EXON;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_INTRON;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.SKIPPED_EXONS;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_1;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.novel.RetainedIntron;
import com.hartwig.hmftools.isofox.novel.RetainedIntronFinder;

import org.junit.Test;

public class NovelJunctionsTest
{
    @Test
    public void testAltSpliceJunctionTypes()
    {
        IsofoxConfig config = createIsofoxConfig();
        config.FragmentSizeData.add(new FragmentSize(30, 1));
        config.ReadLength = 10;

        String chromosome = CHR_1;
        String geneId = GENE_ID_1;

        GeneData geneData = new GeneData(geneId, geneId, chromosome, (byte) 1, 100, 1500, "");

        int transId1 = 1;
        String transName1 = "TRANS01";

        TranscriptData transData1 = new TranscriptData(transId1, transName1, geneId, true, (byte) 1,
                100, 300, null, null, "");

        transData1.exons().add(new ExonData(transId1, 100, 300, 1, -1, -1));
        transData1.exons().add(new ExonData(transId1, 400, 500, 2, -1, -1));
        transData1.exons().add(new ExonData(transId1, 800, 1000, 3, -1, -1));
        transData1.exons().add(new ExonData(transId1, 1400, 1500, 4, -1, -1));

        int transId2 = 2;
        String transName2 = "TRANS02";

        TranscriptData transData2 = new TranscriptData(transId2, transName2, geneId, true, (byte) 1,
                100, 300, null, null, "");

        transData2.exons().add(new ExonData(transId2, 250, 350, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 400, 500, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 600, 700, 3, -1, -1));

        // single exons cannot match a splice junction
        int transId3 = 3;
        String transName3 = "TRANS03";

        TranscriptData transData3 = new TranscriptData(transId3, transName3, geneId, true, (byte) 1,
                100, 300, null, null, "");

        transData3.exons().add(new ExonData(transId3, 100, 300, 1, -1, -1));

        int transId4 = 4;
        String transName4 = "TRANS04";

        TranscriptData transData4 = new TranscriptData(transId4, transName4, geneId, true, (byte) 1,
                360, 500, null, null, "");

        transData4.exons().add(new ExonData(transId4, 360, 500, 1, -1, -1));

        GeneReadData gene = new GeneReadData(geneData);

        List<TranscriptData> transcripts = Lists.newArrayList(transData1, transData2, transData3, transData4);
        List<Integer> transIds = Lists.newArrayList(transId1, transId2, transId3, transId4);

        gene.setTranscripts(transcripts);

        AltSpliceJunctionFinder asjFinder = new AltSpliceJunctionFinder(createIsofoxConfig(), null);
        GeneCollection genes = new GeneCollection(0, Lists.newArrayList(gene));
        asjFinder.setGeneData(genes);

        int[] spliceJunction = new int[SE_PAIR];

        spliceJunction[SE_START] = 300;
        spliceJunction[SE_END] = 400;
        assertTrue(asjFinder.junctionMatchesKnownSpliceJunction(genes.genes(), spliceJunction, transIds));

        spliceJunction[SE_START] = 300;
        spliceJunction[SE_END] = 401;
        assertFalse(asjFinder.junctionMatchesKnownSpliceJunction(genes.genes(), spliceJunction, transIds));

        spliceJunction[SE_START] = 350;
        spliceJunction[SE_END] = 800;
        assertFalse(asjFinder.junctionMatchesKnownSpliceJunction(genes.genes(), spliceJunction, transIds));

        // test out various types of novel junctions:

        // known 5' to novel intronic 3'
        ReadRecord read = createReadRecord(1, chromosome, 291, 369, REF_BASE_STR_1, createCigar(0, 10, 59, 10, 0));

        List<RegionReadData> overlappingRegions = gene.findOverlappingRegions(read);

        read.processOverlappingRegions(overlappingRegions);
        AltSpliceJunction altSJ = asjFinder.createFromRead(read, transIds);

        assertEquals(NOVEL_3_PRIME, altSJ.type());
        assertEquals(300, altSJ.SpliceJunction[SE_START]);
        assertEquals(360, altSJ.SpliceJunction[SE_END]);
        assertEquals(SPLICE_JUNC, altSJ.RegionContexts[SE_START]);
        assertEquals(AltSpliceJunctionContext.INTRONIC, altSJ.RegionContexts[SE_END]);
        assertTrue(altSJ.getSjStartRegions().stream().anyMatch(x -> x.hasTransId(transId1)));
        assertTrue(altSJ.getSjEndRegions().isEmpty());

        List<Integer> validTransIds = altSJ.candidateTransIds();
        assertTrue(validTransIds.contains(transId1));
        assertFalse(validTransIds.contains(transId4));

        // known 5' to novel intronic 3'
        read = createReadRecord(1, chromosome, 491, 559, REF_BASE_STR_1, createCigar(0, 10, 49, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);

        read.processOverlappingRegions(overlappingRegions);
        altSJ = asjFinder.createFromRead(read, transIds);
        altSJ.setGeneId(gene.GeneData.GeneId);

        assertEquals(NOVEL_3_PRIME, altSJ.type());
        assertEquals(0, altSJ.calcNearestExonBoundary(SE_START, gene));
        assertEquals(50, altSJ.calcNearestExonBoundary(SE_END, gene));

        // novel intronic 5' to known 3' prime
        read = createReadRecord(1, chromosome, 361, 409, REF_BASE_STR_1, createCigar(0, 10, 29, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);
        read.processOverlappingRegions(overlappingRegions);
        altSJ = asjFinder.createFromRead(read, transIds);
        altSJ.setGeneId(gene.GeneData.GeneId);

        assertEquals(NOVEL_5_PRIME, altSJ.type());

        // intronic in both
        read = createReadRecord(1, chromosome, 721, 779, REF_BASE_STR_1, createCigar(0, 10, 39, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);
        read.processOverlappingRegions(overlappingRegions);

        assertTrue(overlappingRegions.isEmpty());
        altSJ = asjFinder.createFromRead(read, transIds);
        altSJ.setGeneId(gene.GeneData.GeneId);

        assertEquals(INTRONIC, altSJ.type());
        assertEquals(30, altSJ.calcNearestExonBoundary(SE_START, gene));
        assertEquals(30, altSJ.calcNearestExonBoundary(SE_END, gene));

        // skipped exons
        read = createReadRecord(1, chromosome, 291, 809, REF_BASE_STR_1, createCigar(0, 10, 499, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);
        read.processOverlappingRegions(overlappingRegions);
        altSJ = asjFinder.createFromRead(read, transIds);
        altSJ.setGeneId(gene.GeneData.GeneId);

        assertEquals(SKIPPED_EXONS, altSJ.type());
        assertEquals(0, altSJ.calcNearestExonBoundary(SE_START, gene));
        assertEquals(0, altSJ.calcNearestExonBoundary(SE_END, gene));

        // novel intron
        read = createReadRecord(1, chromosome, 841, 959, REF_BASE_STR_1, createCigar(0, 10, 99, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);
        read.processOverlappingRegions(overlappingRegions);
        altSJ = asjFinder.createFromRead(read, transIds);
        altSJ.setGeneId(gene.GeneData.GeneId);

        assertEquals(NOVEL_INTRON, altSJ.type());
        assertEquals(-150, altSJ.calcNearestExonBoundary(SE_START, gene));
        assertEquals(-150, altSJ.calcNearestExonBoundary(SE_END, gene));

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

        // SJs matching different transcripts
        read = createReadRecord(1, chromosome, 291, 609, REF_BASE_STR_1, createCigar(0, 10, 299, 10, 0));

        overlappingRegions = gene.findOverlappingRegions(read);
        read.processOverlappingRegions(overlappingRegions);
        transIds = read.getTranscriptClassifications().keySet().stream().collect(Collectors.toList());

        altSJ = asjFinder.createFromRead(read, transIds);
        altSJ.setGeneId(gene.GeneData.GeneId);

        assertEquals(MIXED_TRANS, altSJ.type());
        assertEquals(300, altSJ.SpliceJunction[SE_START]);
        assertEquals(600, altSJ.SpliceJunction[SE_END]);
        assertEquals(SPLICE_JUNC, altSJ.RegionContexts[SE_START]);
        assertEquals(SPLICE_JUNC, altSJ.RegionContexts[SE_END]);

        validTransIds = altSJ.candidateTransIds();
        assertTrue(validTransIds.contains(transId1));
        assertTrue(validTransIds.contains(transId2));

        // circular exon looking like a DP
        ReadRecord[] readPair = createSupplementaryReadPair(1, genes, genes, 400, 419, 481, 500,
                createCigar(5, 20, 0), createCigar(0, 20, 5), true);
        readPair[0].processOverlappingRegions(gene.findOverlappingRegions(readPair[0]));
        readPair[1].processOverlappingRegions(gene.findOverlappingRegions(readPair[1]));

        transIds = Lists.newArrayList(transId1);
        AltSpliceJunction circularAltSJ = asjFinder.createFromReads(readPair[0], readPair[1], transIds);

        assertTrue(circularAltSJ != null);
        assertEquals(CIRCULAR, circularAltSJ.type());
        assertEquals(400, circularAltSJ.SpliceJunction[SE_START]);
        assertEquals(500, circularAltSJ.SpliceJunction[SE_END]);
        assertEquals(SPLICE_JUNC, circularAltSJ.RegionContexts[SE_START]);
        assertEquals(SPLICE_JUNC, circularAltSJ.RegionContexts[SE_END]);
    }

    @Test
    public void testRetainedIntrons()
    {
        IsofoxConfig config = createIsofoxConfig();
        config.FragmentSizeData.add(new FragmentSize(30, 1));
        config.ReadLength = 10;

        String chromosome = CHR_1;
        String geneId = GENE_ID_1;

        GeneData geneData = new GeneData(geneId, geneId, chromosome, (byte) 1, 100, 1500, "");

        int transId1 = 1;
        String transName1 = "TRANS01";

        TranscriptData transData1 = new TranscriptData(transId1, transName1, geneId, true, (byte) 1,
                100, 300, null, null, "");

        transData1.exons().add(new ExonData(transId1, 100, 300, 1, -1, -1));
        transData1.exons().add(new ExonData(transId1, 400, 500, 2, -1, -1));
        transData1.exons().add(new ExonData(transId1, 800, 1000, 3, -1, -1));
        transData1.exons().add(new ExonData(transId1, 1400, 1500, 4, -1, -1));

        int transId2 = 2;
        String transName2 = "TRANS02";

        TranscriptData transData2 = new TranscriptData(transId2, transName2, geneId, true, (byte) 1,
                100, 300, null, null, "");

        transData2.exons().add(new ExonData(transId2, 100, 310, 1, -1, -1));
        transData2.exons().add(new ExonData(transId2, 400, 500, 2, -1, -1));
        transData2.exons().add(new ExonData(transId2, 600, 700, 3, -1, -1));

        // single exons cannot match a splice junction
        int transId3 = 3;
        String transName3 = "TRANS03";

        TranscriptData transData3 = new TranscriptData(transId3, transName3, geneId, true, (byte) 1,
                100, 310, null, null, "");

        transData3.exons().add(new ExonData(transId3, 100, 300, 1, -1, -1));
        transData3.exons().add(new ExonData(transId3, 420, 500, 2, -1, -1));

        GeneReadData gene = new GeneReadData(geneData);

        List<TranscriptData> transcripts = Lists.newArrayList(transData1, transData2, transData3);
        List<Integer> transIds = Lists.newArrayList(transId1, transId2, transId3);

        gene.setTranscripts(transcripts);

        RetainedIntronFinder riFinder = new RetainedIntronFinder(config, null);

        GeneCollection genes = new GeneCollection(0, Lists.newArrayList(gene));
        riFinder.setGeneData(genes);

        // first read doesn't span an exon-intron boundary for every transcript
        ReadRecord read1 = createReadRecord(1, chromosome, 291, 310, REF_BASE_STR_1, createCigar(0, 20, 0));
        ReadRecord read2 = createReadRecord(1, chromosome, 340, 360, REF_BASE_STR_1, createCigar(0, 20, 0));

        read1.processOverlappingRegions(gene.findOverlappingRegions(read1));
        read1.processOverlappingRegions(gene.findOverlappingRegions(read2));

        riFinder.evaluateFragmentReads(read1, read2);

        assertEquals(0, riFinder.getRetainedIntrons().size());

        // first read supports a retained intron
        read1 = createReadRecord(1, chromosome, 281, 320, REF_BASE_STR_2, createCigar(0, 40, 0));
        read2 = createReadRecord(1, chromosome, 340, 360, REF_BASE_STR_1, createCigar(0, 20, 0));

        read1.processOverlappingRegions(gene.findOverlappingRegions(read1));
        read1.processOverlappingRegions(gene.findOverlappingRegions(read2));

        riFinder.evaluateFragmentReads(read1, read2);

        assertEquals(1, riFinder.getRetainedIntrons().size());
        RetainedIntron retIntron = riFinder.getRetainedIntrons().get(0);
        assertEquals(false, retIntron.isStart());
        assertEquals(310, retIntron.position());
        assertEquals(1, retIntron.regions().size());
        assertTrue(retIntron.regions().stream().anyMatch(x -> x.hasTransId(transId2)));

        read1 = createReadRecord(1, chromosome, 391, 430, REF_BASE_STR_2, createCigar(0, 40, 0));
        read2 = createReadRecord(1, chromosome, 440, 460, REF_BASE_STR_1, createCigar(0, 20, 0));

        read1.processOverlappingRegions(gene.findOverlappingRegions(read1));
        read1.processOverlappingRegions(gene.findOverlappingRegions(read2));

        riFinder.evaluateFragmentReads(read1, read2);

        assertEquals(2, riFinder.getRetainedIntrons().size());
        retIntron = riFinder.getRetainedIntrons().get(1);
        assertEquals(true, retIntron.isStart());
        assertEquals(400, retIntron.position());
        assertEquals(1, retIntron.regions().size());
        assertTrue(retIntron.regions().stream().anyMatch(x -> x.hasTransId(transId1)));
        assertTrue(retIntron.regions().stream().anyMatch(x -> x.hasTransId(transId2)));

        // with splice support on the other read - updating the same retained intron
        read1 = createReadRecord(1, chromosome, 391, 430, REF_BASE_STR_2, createCigar(0, 40, 0));
        read2 = createReadRecord(1, chromosome, 491, 609, REF_BASE_STR_1, createCigar(0, 10, 99, 10, 0));

        read1.processOverlappingRegions(gene.findOverlappingRegions(read1));
        read2.processOverlappingRegions(gene.findOverlappingRegions(read2));

        riFinder.evaluateFragmentReads(read1, read2);

        assertEquals(2, riFinder.getRetainedIntrons().size());
        retIntron = riFinder.getRetainedIntrons().get(1);
        assertEquals(true, retIntron.isStart());
        assertEquals(400, retIntron.position());
        assertEquals(1, retIntron.regions().size());
        assertEquals(2, retIntron.getFragmentCount());
        assertEquals(1, retIntron.getSplicedFragmentCount());
        assertTrue(retIntron.regions().stream().anyMatch(x -> x.hasTransId(transId1)));
        assertTrue(retIntron.regions().stream().anyMatch(x -> x.hasTransId(transId2)));

        // cannot be the first or last exon of a transcript
        read1 = createReadRecord(1, chromosome, 91, 110, REF_BASE_STR_1, createCigar(0, 20, 0));
        read2 = createReadRecord(1, chromosome, 121, 140, REF_BASE_STR_1, createCigar(0, 20, 0));

        read1.processOverlappingRegions(gene.findOverlappingRegions(read1));
        read2.processOverlappingRegions(gene.findOverlappingRegions(read2));

        riFinder.evaluateFragmentReads(read1, read2);

        assertEquals(2, riFinder.getRetainedIntrons().size());

        read1 = createReadRecord(1, chromosome, 1491, 1510, REF_BASE_STR_1, createCigar(0, 20, 0));
        read2 = createReadRecord(1, chromosome, 1551, 1570, REF_BASE_STR_1, createCigar(0, 20, 0));

        read1.processOverlappingRegions(gene.findOverlappingRegions(read1));
        read2.processOverlappingRegions(gene.findOverlappingRegions(read2));

        riFinder.evaluateFragmentReads(read1, read2);

        assertEquals(2, riFinder.getRetainedIntrons().size());

        riFinder.getRetainedIntrons().clear();

        // both reads covering boundaries of the same exon - not permitted
        read1 = createReadRecord(1, chromosome, 391, 410, REF_BASE_STR_1, createCigar(0, 20, 0));
        read2 = createReadRecord(1, chromosome, 491, 510, REF_BASE_STR_1, createCigar(0, 20, 0));

        read1.processOverlappingRegions(gene.findOverlappingRegions(read1));
        read2.processOverlappingRegions(gene.findOverlappingRegions(read2));

        riFinder.evaluateFragmentReads(read1, read2);

        assertTrue(riFinder.getRetainedIntrons().isEmpty());


        // both reads covering boundaries of sequential exons is permitted
        read1 = createReadRecord(1, chromosome, 491, 510, REF_BASE_STR_1, createCigar(0, 20, 0));
        read2 = createReadRecord(1, chromosome, 591, 610, REF_BASE_STR_1, createCigar(0, 20, 0));

        read1.processOverlappingRegions(gene.findOverlappingRegions(read1));
            read2.processOverlappingRegions(gene.findOverlappingRegions(read2));

        riFinder.evaluateFragmentReads(read1, read2);

        assertEquals(2, riFinder.getRetainedIntrons().size());

        retIntron = riFinder.getRetainedIntrons().get(0);
        assertEquals(false, retIntron.isStart());
        assertEquals(500, retIntron.position());
        assertEquals(1, retIntron.regions().size());
        assertEquals(1, retIntron.getFragmentCount());
        assertTrue(retIntron.regions().stream().anyMatch(x -> x.hasTransId(transId1)));
        assertTrue(retIntron.regions().stream().anyMatch(x -> x.hasTransId(transId2)));

        retIntron = riFinder.getRetainedIntrons().get(1);
        assertEquals(true, retIntron.isStart());
        assertEquals(600, retIntron.position());
        assertEquals(1, retIntron.regions().size());
        assertEquals(1, retIntron.getFragmentCount());
    }

}
