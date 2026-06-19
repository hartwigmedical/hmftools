package com.hartwig.hmftools.tars.liftback.rescue;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.junit.Test;

public class AnnotatedJunctionLoaderTest
{
    private static final byte POS_STRAND = 1;
    private static final byte NEG_STRAND = -1;
    private static final String BIOTYPE = "protein_coding";

    private static GeneData gene(final String geneId, final String chromosome)
    {
        return createEnsemblGeneData(geneId, geneId, chromosome, 1, 1, 1_000_000);
    }

    // exonLength 99 -> a 100-base exon spanning start..start+99 (e.g. 1000..1099).
    private static TranscriptData transcript(final String geneId, final int transId, final byte strand, final int[] exonStarts)
    {
        return transcript(geneId, transId, strand, exonStarts, 99);
    }

    private static TranscriptData transcript(
            final String geneId, final int transId, final byte strand, final int[] exonStarts, final int exonLength)
    {
        return createTransExons(geneId, transId, strand, exonStarts, exonLength, null, null, false, BIOTYPE);
    }

    @Test
    public void testLoadsIntronsFromAdjacentExonPairs()
    {
        final EnsemblDataCache cache = createGeneDataCache();
        addGeneData(cache, "1", List.of(gene("ENSG_TEST", "1")));

        // TX_A: 3 exons -> 2 junctions; TX_B: different middle exon -> 1 distinct junction
        addTransExonData(cache, "ENSG_TEST", List.of(
                transcript("ENSG_TEST", 1, POS_STRAND, new int[] { 1000, 1500, 2000 }),
                transcript("ENSG_TEST", 2, POS_STRAND, new int[] { 1000, 1700 })));

        final Set<ChrBaseRegion> introns = AnnotatedJunctionLoader.deriveIntrons(cache, V38);

        assertEquals(3, introns.size());
        assertTrue(introns.contains(new ChrBaseRegion("chr1", 1100, 1499)));
        assertTrue(introns.contains(new ChrBaseRegion("chr1", 1600, 1999)));
        assertTrue(introns.contains(new ChrBaseRegion("chr1", 1100, 1699)));
    }

    @Test
    public void testChrPrefixNormalization()
    {
        final EnsemblDataCache cache = createGeneDataCache();
        addGeneData(cache, "chr3", List.of(gene("ENSG_PREFIXED", "chr3")));
        addGeneData(cache, "3", List.of(gene("ENSG_BARE", "3")));

        addTransExonData(cache, "ENSG_PREFIXED", List.of(transcript("ENSG_PREFIXED", 1, POS_STRAND, new int[] { 100, 500 })));
        addTransExonData(cache, "ENSG_BARE", List.of(transcript("ENSG_BARE", 2, POS_STRAND, new int[] { 1000, 1500 })));

        final Set<ChrBaseRegion> introns = AnnotatedJunctionLoader.deriveIntrons(cache, V38);

        // both genes' introns land on chr3 regardless of the cache's chr-prefix form
        assertEquals(2, introns.size());
        assertTrue(introns.contains(new ChrBaseRegion("chr3", 200, 499)));
        assertTrue(introns.contains(new ChrBaseRegion("chr3", 1100, 1499)));
    }

    @Test
    public void testSingleExonTranscriptsProduceNoIntrons()
    {
        final EnsemblDataCache cache = createGeneDataCache();
        addGeneData(cache, "1", List.of(gene("ENSG_X", "1")));
        addTransExonData(cache, "ENSG_X", List.of(transcript("ENSG_X", 1, POS_STRAND, new int[] { 1000 }, 1000)));

        final Set<ChrBaseRegion> introns = AnnotatedJunctionLoader.deriveIntrons(cache, V38);

        assertEquals(0, introns.size());
    }

    @Test
    public void testDuplicateJunctionsAcrossTranscriptsDeduplicate()
    {
        final EnsemblDataCache cache = createGeneDataCache();
        addGeneData(cache, "2", List.of(gene("ENSG_DUP", "2")));
        addTransExonData(cache, "ENSG_DUP", List.of(
                transcript("ENSG_DUP", 1, POS_STRAND, new int[] { 100, 300 }),
                transcript("ENSG_DUP", 2, POS_STRAND, new int[] { 100, 300 })));

        final Set<ChrBaseRegion> introns = AnnotatedJunctionLoader.deriveIntrons(cache, V38);

        assertEquals(1, introns.size());
        assertTrue(introns.contains(new ChrBaseRegion("chr2", 200, 299)));
    }

    @Test
    public void testReverseStrandRankOrder()
    {
        // negative-strand transcripts have exons ranked high->low genome pos; intron coords must still be lo..hi.
        final EnsemblDataCache cache = createGeneDataCache();
        addGeneData(cache, "4", List.of(gene("ENSG_NEG", "4")));
        addTransExonData(cache, "ENSG_NEG", List.of(transcript("ENSG_NEG", 1, NEG_STRAND, new int[] { 1000, 2000 })));

        final Set<ChrBaseRegion> introns = AnnotatedJunctionLoader.deriveIntrons(cache, V38);

        assertEquals(1, introns.size());
        assertTrue(introns.contains(new ChrBaseRegion("chr4", 1100, 1999)));
    }
}
