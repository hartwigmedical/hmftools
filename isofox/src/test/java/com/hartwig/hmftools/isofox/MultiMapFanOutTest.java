package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.isofox.FragmentAllocator.MULTI_MAP_SPLICED;
import static com.hartwig.hmftools.isofox.FragmentAllocator.MULTI_MAP_UNSPLICED;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.TestUtils.ALT_SJ_COHORT_CACHE;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_3;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_4;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createGeneCollection;
import static com.hartwig.hmftools.isofox.TestUtils.createIsofoxConfig;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.fusion.FusionTestUtils.createGeneDataCache;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.Read;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

import org.junit.Test;

// exercises the multi-map fan-out numerics (exon gate, mass renormalisation, spliced bucketing, divergent-mate
// pinning) end-to-end through processReadRecords; MultiMapReadTest separately covers XA-string parsing
public class MultiMapFanOutTest
{
    private static final double EPSILON = 1e-6;

    // GENE4 exon1 on chr2 (1000-1100); an alt landing here is exon-supported
    private static final ChrBaseRegion GENE4_EXON = new ChrBaseRegion(CHR_2, 1000, 1099);

    // GENE3 exon1 on chr1 (20000-20100); another exon-supported locus
    private static final ChrBaseRegion GENE3_EXON = new ChrBaseRegion(CHR_1, 20000, 20099);

    // between GENE2 exon1 (10000-10100) and exon2 (10200-10300); overlaps the gene span but no exon
    private static final ChrBaseRegion GENE2_INTRON = new ChrBaseRegion(CHR_1, 10120, 10180);

    private static FragmentAllocator createAllocator(final EnsemblDataCache geneTransCache)
    {
        IsofoxConfig config = createIsofoxConfig();
        config.Functions.clear();
        config.Functions.add(TRANSCRIPT_COUNTS);
        return new FragmentAllocator(config, geneTransCache, ALT_SJ_COHORT_CACHE, new ResultsWriter(config));
    }

    // a concordant, non-chimeric pair mapped inside GENE1 on chr1; the caller then attaches XA alt loci
    private static Read[] createGene1Pair()
    {
        Read read1 = createReadRecord(1, CHR_1, 1000, 1099, generateRandomBases(100), createCigar(0, 100, 0));
        read1.setFragmentInsertSize(500);

        Read read2 = createReadRecord(1, CHR_1, 1400, 1499, generateRandomBases(100), createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-500);

        return new Read[] { read1, read2 };
    }

    private static List<Read.AltAlignment> altLoci(final Read.AltAlignment... alts)
    {
        return Lists.newArrayList(alts);
    }

    @Test
    public void testExonGateDropsIntronicAltAndRenormalises()
    {
        EnsemblDataCache geneTransCache = createGeneDataCache();
        FragmentAllocator allocator = createAllocator(geneTransCache);
        GeneCollection geneSet = createGeneCollection(geneTransCache, 0, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));

        // one exon-supported alt (GENE4) and one intronic alt (GENE2); both mates carry the same two alts
        Read[] reads = createGene1Pair();
        List<Read.AltAlignment> alts = altLoci(
                new Read.AltAlignment(GENE4_EXON, false),
                new Read.AltAlignment(GENE2_INTRON, false));
        reads[0].setAltLoci(alts);
        reads[1].setAltLoci(alts);

        allocator.processReadRecords(geneSet, Lists.newArrayList(reads));

        Map<String,double[]> counts = allocator.getMultiMapGeneCounts();

        // the intronic alt supports no transcript so it is not credited and does not inflate the locus count;
        // the fragment renormalises over the surviving loci (primary + GENE4) so GENE4 gets 1/2, not 1/3
        assertFalse(counts.containsKey(GENE_ID_2));
        assertTrue(counts.containsKey(GENE_ID_4));
        assertEquals(0.5, counts.get(GENE_ID_4)[MULTI_MAP_UNSPLICED], EPSILON);
        assertEquals(0.0, counts.get(GENE_ID_4)[MULTI_MAP_SPLICED], EPSILON);
    }

    @Test
    public void testTwoExonicAltsSplitOverThreeLoci()
    {
        EnsemblDataCache geneTransCache = createGeneDataCache();
        FragmentAllocator allocator = createAllocator(geneTransCache);
        GeneCollection geneSet = createGeneCollection(geneTransCache, 0, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));

        // two exon-supported alts: GENE4 (unspliced) and GENE3 (spliced); denominator = primary + 2 alts = 3
        Read[] reads = createGene1Pair();
        List<Read.AltAlignment> alts = altLoci(
                new Read.AltAlignment(GENE4_EXON, false),
                new Read.AltAlignment(GENE3_EXON, true));
        reads[0].setAltLoci(alts);
        reads[1].setAltLoci(alts);

        allocator.processReadRecords(geneSet, Lists.newArrayList(reads));

        Map<String,double[]> counts = allocator.getMultiMapGeneCounts();

        assertEquals(1.0 / 3, counts.get(GENE_ID_4)[MULTI_MAP_UNSPLICED], EPSILON);
        assertEquals(0.0, counts.get(GENE_ID_4)[MULTI_MAP_SPLICED], EPSILON);

        // GENE3's alt CIGAR was spliced, so its share is booked to the spliced bucket
        assertEquals(1.0 / 3, counts.get(GENE_ID_3)[MULTI_MAP_SPLICED], EPSILON);
        assertEquals(0.0, counts.get(GENE_ID_3)[MULTI_MAP_UNSPLICED], EPSILON);
    }

    @Test
    public void testUniqueMatePinsFragment()
    {
        EnsemblDataCache geneTransCache = createGeneDataCache();
        FragmentAllocator allocator = createAllocator(geneTransCache);
        GeneCollection geneSet = createGeneCollection(geneTransCache, 0, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));

        // read1 is multi-mapped but read2 is uniquely placed; the unique mate pins the fragment, so no fan-out
        Read[] reads = createGene1Pair();
        reads[0].setAltLoci(altLoci(new Read.AltAlignment(GENE4_EXON, false)));

        allocator.processReadRecords(geneSet, Lists.newArrayList(reads));

        assertTrue(allocator.getMultiMapGeneCounts().isEmpty());
    }

    @Test
    public void testAllIntronicAltsCreditNothing()
    {
        EnsemblDataCache geneTransCache = createGeneDataCache();
        FragmentAllocator allocator = createAllocator(geneTransCache);
        GeneCollection geneSet = createGeneCollection(geneTransCache, 0, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));

        // both alts land in introns; nothing is credited and the mass stays on the primary locus
        Read[] reads = createGene1Pair();
        List<Read.AltAlignment> alts = altLoci(
                new Read.AltAlignment(GENE2_INTRON, false),
                new Read.AltAlignment(new ChrBaseRegion(CHR_1, 10320, 10380), false)); // GENE2 intron 2/3
        reads[0].setAltLoci(alts);
        reads[1].setAltLoci(alts);

        allocator.processReadRecords(geneSet, Lists.newArrayList(reads));

        assertTrue(allocator.getMultiMapGeneCounts().isEmpty());
    }
}
