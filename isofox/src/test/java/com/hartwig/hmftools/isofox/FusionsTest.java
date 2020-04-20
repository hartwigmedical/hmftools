package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.generateExonStarts;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_1;
import static com.hartwig.hmftools.isofox.ReadCountsTest.createCigar;
import static com.hartwig.hmftools.isofox.ReadCountsTest.createReadRecord;
import static com.hartwig.hmftools.isofox.common.ReadRecord.findOverlappingRegions;
import static com.hartwig.hmftools.isofox.common.RegionReadData.generateExonicRegions;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.SPLICED_BOTH;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.fusion.FusionFragment;

import org.junit.Test;

import htsjdk.samtools.Cigar;

public class FusionsTest
{
    private static final String CHR_1 = "1";

    private static final String GENE_NAME_1 = "GENE1"; // +ve strand
    private static final String GENE_ID_1 = "ENSG0001";
    private static final long GENE_START_1 = 1000;

    private static final String GENE_NAME_2 = "GENE2"; // +ve strand
    private static final String GENE_ID_2 = "ENSG0002";
    private static final long GENE_START_2 = 10000;

    private static final String GENE_NAME_3 = "GENE3"; // -ve strand
    private static final String GENE_ID_3 = "ENSG0003";
    private static final long GENE_START_3 = 20000;

    private static final String CHR_2 = "2";

    private static final String GENE_NAME_4 = "GENE4"; // +ve strand
    private static final String GENE_ID_4 = "ENSG0004";
    private static final long GENE_START_4 = 1000;

    private static final String GENE_NAME_5 = "GENE5"; // -ve strand
    private static final String GENE_ID_5 = "ENSG0005";
    private static final long GENE_START_5 = 10000;

    private static final byte POS_STRAND = 1;
    private static final byte NEG_STRAND = -1;

    private static final int EXON_LENGTH = 100;
    private static final int INTRON_LENGTH = 200;

    private static int getGeneCollection(final String geneId)
    {
        if(geneId == GENE_ID_1)
            return 0;
        if(geneId == GENE_ID_2)
            return 1;
        if(geneId == GENE_ID_3)
            return 2;
        if(geneId == GENE_ID_4)
            return 3;
        if(geneId == GENE_ID_5)
            return 4;

        return -1;
    }

    private static void addTestGenes(EnsemblDataCache geneTransCache)
    {
        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, GENE_START_1, GENE_START_1 + 1000));
        geneList.add(createEnsemblGeneData(GENE_ID_2, GENE_NAME_2, CHR_1, POS_STRAND, GENE_START_2, GENE_START_2 + 1000));
        geneList.add(createEnsemblGeneData(GENE_ID_3, GENE_NAME_3, CHR_1, NEG_STRAND, GENE_START_3, GENE_START_3 + 1000));
        addGeneData(geneTransCache, CHR_1, geneList);

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(GENE_ID_4, GENE_NAME_4, CHR_2, POS_STRAND, GENE_START_4, GENE_START_4 + 1000));
        geneList.add(createEnsemblGeneData(GENE_ID_5, GENE_NAME_5, CHR_2, NEG_STRAND, GENE_START_5, GENE_START_5 + 1000));
        addGeneData(geneTransCache, CHR_2, geneList);
    }

    private static final int TRANS_1 = 1;
    private static final int TRANS_2 = 2;
    private static final int TRANS_3 = 3;
    private static final int TRANS_4 = 4;
    private static final int TRANS_5 = 5;

    private static void addTestTranscripts(EnsemblDataCache geneTransCache)
    {
        Long codingStart = null;
        Long codingEnd = null;
        boolean canonical = true;

        List<TranscriptData> transDataList = Lists.newArrayList();

        TranscriptData transData = createTransExons(
                GENE_ID_1, TRANS_1, POS_STRAND, generateExonStarts(GENE_START_1, 3, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_1, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_2, TRANS_2, POS_STRAND, generateExonStarts(GENE_START_2, 3, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_2, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_3, TRANS_3, NEG_STRAND, generateExonStarts(GENE_START_3, 3, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_3, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_4, TRANS_4, POS_STRAND, generateExonStarts(GENE_START_4, 3, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_4, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_5, TRANS_5, NEG_STRAND, generateExonStarts(GENE_START_5, 3, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_5, transDataList);
    }

    private static ReadRecord createMappedRead(final int id, final GeneCollection geneCollection, long posStart, long posEnd, final Cigar cigar)
    {
        ReadRecord read = createReadRecord(id, geneCollection.chromosome(), posStart, posEnd, REF_BASE_STR_1, cigar);

        read.processOverlappingRegions(findOverlappingRegions(geneCollection.getExonRegions(), read));
        read.captureGeneInfo(geneCollection.id());

        return read;
    }

    private GeneCollection createGeneCollection(final EnsemblDataCache geneTransCache, int id, final List<EnsemblGeneData> genes)
    {
        List<GeneReadData> geneReadData = Lists.newArrayList();
        genes.forEach(x -> geneReadData.add(new GeneReadData(x)));
        geneReadData.forEach(x -> x.setTranscripts(geneTransCache.getTranscripts(x.GeneData.GeneId)));

        GeneCollection gc = new GeneCollection(id, geneReadData);
        return gc;
    }

    @Test
    public void testFragmentClassification()
    {
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        // a simple DEL
        ReadRecord read1 = createMappedRead(1, gc1, 1050, 1089, createCigar(0, 40, 0));
        ReadRecord read2 = createMappedRead(1, gc1, 1081, 1100, createCigar(0, 20, 20));
        ReadRecord read3 = createMappedRead(1, gc2, 10200, 10219, createCigar(20, 20, 0));

        List<ReadRecord> reads = Lists.newArrayList(read1, read2, read3);
        FusionFragment fragment = new FusionFragment(reads);

        assertEquals(SPLICED_BOTH, fragment.type());
        assertEquals(CHR_1, fragment.chromosomes()[SE_START]);
        assertEquals(CHR_1, fragment.chromosomes()[SE_END]);
        assertEquals(1100, fragment.splicePositions()[SE_START]);
        assertEquals(10200, fragment.splicePositions()[SE_END]);
        assertEquals(1, fragment.spliceOrientations()[SE_START]);
        assertEquals(-1, fragment.spliceOrientations()[SE_END]);
        assertEquals(SPLICED_BOTH, fragment.type());
        assertTrue(fragment.hasValidSpliceData());
        assertEquals(DEL, fragment.getImpliedSvType());

        final List<List<String>> spliceGeneIds = Lists.newArrayList(Lists.newArrayList(), Lists.newArrayList());
        fragment.populateGeneCandidates(spliceGeneIds);

        assertEquals(1, spliceGeneIds.get(SE_START).size());
        assertEquals(GENE_ID_1, spliceGeneIds.get(SE_START).get(0));
        assertEquals(1, spliceGeneIds.get(SE_END).size());
        assertEquals(GENE_ID_2, spliceGeneIds.get(SE_END).get(0));

        // DUP
        read1 = createMappedRead(1, gc2, 10220, 10259, createCigar(0, 40, 0));
        read2 = createMappedRead(1, gc2, 10281, 10300, createCigar(0, 20, 20));
        read3 = createMappedRead(1, gc1, 1200, 1219, createCigar(20, 20, 0));

        reads = Lists.newArrayList(read1, read2, read3);
        fragment = new FusionFragment(reads);

        assertEquals(SPLICED_BOTH, fragment.type());
        assertEquals(CHR_1, fragment.chromosomes()[SE_START]);
        assertEquals(CHR_1, fragment.chromosomes()[SE_END]);
        assertEquals(1200, fragment.splicePositions()[SE_START]);
        assertEquals(10300, fragment.splicePositions()[SE_END]);
        assertEquals(-1, fragment.spliceOrientations()[SE_START]);
        assertEquals(1, fragment.spliceOrientations()[SE_END]);
        assertEquals(SPLICED_BOTH, fragment.type());
        assertTrue(fragment.hasValidSpliceData());
        assertEquals(DUP, fragment.getImpliedSvType());

        fragment.populateGeneCandidates(spliceGeneIds);

        assertEquals(1, spliceGeneIds.get(SE_START).size());
        assertEquals(GENE_ID_1, spliceGeneIds.get(SE_START).get(0));
        assertEquals(1, spliceGeneIds.get(SE_END).size());
        assertEquals(GENE_ID_2, spliceGeneIds.get(SE_END).get(0));

        // INV
        final GeneCollection gc3 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_3)));

        read1 = createMappedRead(1, gc1, 1081, 1100, createCigar(0, 20, 20));
        read2 = createMappedRead(1, gc3, 20281, 20300, createCigar(0, 20, 20));
        read3 = createMappedRead(1, gc3, 20220, 20259, createCigar(20, 20, 0));

        reads = Lists.newArrayList(read1, read2, read3);
        fragment = new FusionFragment(reads);

        assertEquals(SPLICED_BOTH, fragment.type());
        assertEquals(CHR_1, fragment.chromosomes()[SE_START]);
        assertEquals(CHR_1, fragment.chromosomes()[SE_END]);
        assertEquals(1100, fragment.splicePositions()[SE_START]);
        assertEquals(20300, fragment.splicePositions()[SE_END]);
        assertEquals(1, fragment.spliceOrientations()[SE_START]);
        assertEquals(1, fragment.spliceOrientations()[SE_END]);
        assertEquals(SPLICED_BOTH, fragment.type());
        assertTrue(fragment.hasValidSpliceData());
        assertEquals(INV, fragment.getImpliedSvType());

        fragment.populateGeneCandidates(spliceGeneIds);

        assertEquals(1, spliceGeneIds.get(SE_START).size());
        assertEquals(GENE_ID_1, spliceGeneIds.get(SE_START).get(0));
        assertEquals(1, spliceGeneIds.get(SE_END).size());
        assertEquals(GENE_ID_3, spliceGeneIds.get(SE_END).get(0));

        // BND
        final GeneCollection gc5 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_5)));

        read1 = createMappedRead(1, gc5, 10281, 10300, createCigar(0, 20, 20));
        read2 = createMappedRead(1, gc3, 20281, 20300, createCigar(0, 20, 20));
        read3 = createMappedRead(1, gc3, 20220, 20259, createCigar(20, 20, 0));

        reads = Lists.newArrayList(read1, read2, read3);
        fragment = new FusionFragment(reads);

        assertEquals(SPLICED_BOTH, fragment.type());
        assertEquals(CHR_1, fragment.chromosomes()[SE_START]);
        assertEquals(CHR_2, fragment.chromosomes()[SE_END]);
        assertEquals(20300, fragment.splicePositions()[SE_START]);
        assertEquals(10300, fragment.splicePositions()[SE_END]);
        assertEquals(1, fragment.spliceOrientations()[SE_START]);
        assertEquals(1, fragment.spliceOrientations()[SE_END]);
        assertEquals(SPLICED_BOTH, fragment.type());
        assertTrue(fragment.hasValidSpliceData());
        assertEquals(BND, fragment.getImpliedSvType());

        fragment.populateGeneCandidates(spliceGeneIds);

        assertEquals(1, spliceGeneIds.get(SE_START).size());
        assertEquals(GENE_ID_3, spliceGeneIds.get(SE_START).get(0));
        assertEquals(1, spliceGeneIds.get(SE_END).size());
        assertEquals(GENE_ID_5, spliceGeneIds.get(SE_END).get(0));
    }
}
