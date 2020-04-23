package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_3;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_5;
import static com.hartwig.hmftools.isofox.TestUtils.createGeneCollection;
import static com.hartwig.hmftools.isofox.TestUtils.addTestGenes;
import static com.hartwig.hmftools.isofox.TestUtils.addTestTranscripts;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createMappedRead;
import static com.hartwig.hmftools.isofox.common.ReadRecord.findOverlappingRegions;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.BOTH_JUNCTIONS;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.fusion.FusionFinder;
import com.hartwig.hmftools.isofox.fusion.FusionFragment;

import org.junit.Test;

import htsjdk.samtools.Cigar;

public class FusionFragmentsTest
{
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

        assertEquals(BOTH_JUNCTIONS, fragment.type());
        assertEquals(CHR_1, fragment.chromosomes()[SE_START]);
        assertEquals(CHR_1, fragment.chromosomes()[SE_END]);
        assertEquals(1100, fragment.splicePositions()[SE_START]);
        assertEquals(10200, fragment.splicePositions()[SE_END]);
        assertEquals(1, fragment.spliceOrientations()[SE_START]);
        assertEquals(-1, fragment.spliceOrientations()[SE_END]);
        assertTrue(fragment.hasValidSpliceData());
        assertEquals(DEL, fragment.getImpliedSvType());

        final List<String>[] spliceGeneIds = new List[] { Lists.newArrayList(), Lists.newArrayList() };

        for(int se = SE_START; se <= SE_END; ++se)
        {
            spliceGeneIds[se].addAll(fragment.getGeneIds(se));
        }

        assertTrue(fragment.isSpliced());
        assertEquals(1, spliceGeneIds[SE_START].size());
        assertEquals(GENE_ID_1, spliceGeneIds[SE_START].get(0));
        assertEquals(1, spliceGeneIds[SE_END].size());
        assertEquals(GENE_ID_2, spliceGeneIds[SE_END].get(0));

        // unspliced DEL - imagining an SV at 1150 to 10150
        read1 = createMappedRead(1, gc1, 1110, 1149, createCigar(0, 40, 0));
        read2 = createMappedRead(1, gc1, 1131, 1150, createCigar(0, 20, 20));
        read3 = createMappedRead(1, gc2, 10150, 10169, createCigar(20, 20, 0));

        reads = Lists.newArrayList(read1, read2, read3);
        read1.addIntronicTranscriptRefs(gc1.getTranscripts());
        read2.addIntronicTranscriptRefs(gc1.getTranscripts());
        read3.addIntronicTranscriptRefs(gc2.getTranscripts());

        fragment = new FusionFragment(reads);

        assertEquals(BOTH_JUNCTIONS, fragment.type());
        assertEquals(CHR_1, fragment.chromosomes()[SE_START]);
        assertEquals(CHR_1, fragment.chromosomes()[SE_END]);
        assertEquals(1150, fragment.splicePositions()[SE_START]);
        assertEquals(10150, fragment.splicePositions()[SE_END]);
        assertEquals(1, fragment.spliceOrientations()[SE_START]);
        assertEquals(-1, fragment.spliceOrientations()[SE_END]);
        assertTrue(fragment.hasValidSpliceData());
        assertEquals(DEL, fragment.getImpliedSvType());

        for(int se = SE_START; se <= SE_END; ++se)
        {
            spliceGeneIds[se].clear();
            spliceGeneIds[se].addAll(fragment.getGeneIds(se));
        }

        assertTrue(fragment.isUnspliced());
        assertEquals(1, spliceGeneIds[SE_START].size());
        assertEquals(GENE_ID_1, spliceGeneIds[SE_START].get(0));
        assertEquals(1, spliceGeneIds[SE_END].size());
        assertEquals(GENE_ID_2, spliceGeneIds[SE_END].get(0));

        // DUP
        read1 = createMappedRead(1, gc2, 10220, 10259, createCigar(0, 40, 0));
        read2 = createMappedRead(1, gc2, 10281, 10300, createCigar(0, 20, 20));
        read3 = createMappedRead(1, gc1, 1200, 1219, createCigar(20, 20, 0));

        reads = Lists.newArrayList(read1, read2, read3);
        fragment = new FusionFragment(reads);

        assertEquals(BOTH_JUNCTIONS, fragment.type());
        assertEquals(CHR_1, fragment.chromosomes()[SE_START]);
        assertEquals(CHR_1, fragment.chromosomes()[SE_END]);
        assertEquals(1200, fragment.splicePositions()[SE_START]);
        assertEquals(10300, fragment.splicePositions()[SE_END]);
        assertEquals(-1, fragment.spliceOrientations()[SE_START]);
        assertEquals(1, fragment.spliceOrientations()[SE_END]);
        assertTrue(fragment.hasValidSpliceData());
        assertEquals(DUP, fragment.getImpliedSvType());

        for(int se = SE_START; se <= SE_END; ++se)
        {
            spliceGeneIds[se].clear();
            spliceGeneIds[se].addAll(fragment.getGeneIds(se));
        }

        assertEquals(1, spliceGeneIds[SE_START].size());
        assertEquals(GENE_ID_1, spliceGeneIds[SE_START].get(0));
        assertEquals(1, spliceGeneIds[SE_END].size());
        assertEquals(GENE_ID_2, spliceGeneIds[SE_END].get(0));

        // INV
        final GeneCollection gc3 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_3)));

        read1 = createMappedRead(1, gc1, 1081, 1100, createCigar(0, 20, 20));
        read2 = createMappedRead(1, gc3, 20281, 20300, createCigar(0, 20, 20));
        read3 = createMappedRead(1, gc3, 20220, 20259, createCigar(20, 20, 0));

        reads = Lists.newArrayList(read1, read2, read3);
        fragment = new FusionFragment(reads);

        assertEquals(BOTH_JUNCTIONS, fragment.type());
        assertEquals(CHR_1, fragment.chromosomes()[SE_START]);
        assertEquals(CHR_1, fragment.chromosomes()[SE_END]);
        assertEquals(1100, fragment.splicePositions()[SE_START]);
        assertEquals(20300, fragment.splicePositions()[SE_END]);
        assertEquals(1, fragment.spliceOrientations()[SE_START]);
        assertEquals(1, fragment.spliceOrientations()[SE_END]);
        assertTrue(fragment.hasValidSpliceData());
        assertEquals(INV, fragment.getImpliedSvType());

        for(int se = SE_START; se <= SE_END; ++se)
        {
            spliceGeneIds[se].clear();
            spliceGeneIds[se].addAll(fragment.getGeneIds(se));
        }

        assertEquals(1, spliceGeneIds[SE_START].size());
        assertEquals(GENE_ID_1, spliceGeneIds[SE_START].get(0));
        assertEquals(1, spliceGeneIds[SE_END].size());
        assertEquals(GENE_ID_3, spliceGeneIds[SE_END].get(0));

        // BND
        final GeneCollection gc5 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_5)));

        read1 = createMappedRead(1, gc5, 10281, 10300, createCigar(0, 20, 20));
        read2 = createMappedRead(1, gc3, 20281, 20300, createCigar(0, 20, 20));
        read3 = createMappedRead(1, gc3, 20220, 20259, createCigar(20, 20, 0));

        reads = Lists.newArrayList(read1, read2, read3);
        fragment = new FusionFragment(reads);

        assertEquals(BOTH_JUNCTIONS, fragment.type());
        assertEquals(CHR_1, fragment.chromosomes()[SE_START]);
        assertEquals(CHR_2, fragment.chromosomes()[SE_END]);
        assertEquals(20300, fragment.splicePositions()[SE_START]);
        assertEquals(10300, fragment.splicePositions()[SE_END]);
        assertEquals(1, fragment.spliceOrientations()[SE_START]);
        assertEquals(1, fragment.spliceOrientations()[SE_END]);
        assertTrue(fragment.hasValidSpliceData());
        assertEquals(BND, fragment.getImpliedSvType());

        for(int se = SE_START; se <= SE_END; ++se)
        {
            spliceGeneIds[se].clear();
            spliceGeneIds[se].addAll(fragment.getGeneIds(se));
        }

        assertEquals(1, spliceGeneIds[SE_START].size());
        assertEquals(GENE_ID_3, spliceGeneIds[SE_START].get(0));
        assertEquals(1, spliceGeneIds[SE_END].size());
        assertEquals(GENE_ID_5, spliceGeneIds[SE_END].get(0));
    }
}
