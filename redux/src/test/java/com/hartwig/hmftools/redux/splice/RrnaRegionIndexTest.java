package com.hartwig.hmftools.redux.splice;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.junit.Test;

// covers RrnaRegionIndex build (biotype filter, merge of overlapping gene ranges) and overlaps lookup boundaries.
public class RrnaRegionIndexTest
{
    private static final byte POS_STRAND = 1;

    @Test
    public void testBuildPicksUpRrnaAndMtRrnaBiotypes()
    {
        final EnsemblDataCache cache = new EnsemblDataCache("", RefGenomeVersion.V38);

        addGene(cache, "ENSG_RRNA", "RRNA_GENE", "1", 100, 200);
        addTranscript(cache, "ENSG_RRNA", "ENST_RRNA", "rRNA");

        addGene(cache, "ENSG_MT", "MT_RRNA_GENE", "MT", 1000, 1500);
        addTranscript(cache, "ENSG_MT", "ENST_MT", "Mt_rRNA");

        addGene(cache, "ENSG_CODING", "CODING_GENE", "1", 300, 400);
        addTranscript(cache, "ENSG_CODING", "ENST_CODING", "protein_coding");

        addGene(cache, "ENSG_PSEUDO", "RRNA_PSEUDO", "1", 500, 600);
        addTranscript(cache, "ENSG_PSEUDO", "ENST_PSEUDO", "rRNA_pseudogene");

        final RrnaRegionIndex index = RrnaRegionIndex.build(cache, RefGenomeVersion.V37);

        assertEquals(2, index.geneCount());
        assertTrue(index.overlaps("1", 100, 100));
        assertTrue(index.overlaps("1", 200, 200));
        assertTrue(index.overlaps("1", 50, 150));
        assertTrue(index.overlaps("1", 150, 250));
        assertTrue(index.overlaps("MT", 1200, 1300));

        // pseudogene biotype must not match — common confusion point
        assertFalse(index.overlaps("1", 500, 600));
        // coding gene unaffected
        assertFalse(index.overlaps("1", 350, 360));
        // outside any rRNA region
        assertFalse(index.overlaps("1", 201, 299));
        assertFalse(index.overlaps("1", 50, 99));
        // unknown chromosome
        assertFalse(index.overlaps("99", 100, 200));
    }

    @Test
    public void testBuildMergesOverlappingAdjacentRrnaGenes()
    {
        final EnsemblDataCache cache = new EnsemblDataCache("", RefGenomeVersion.V38);

        addGene(cache, "ENSG_A", "A", "5", 100, 200);
        addTranscript(cache, "ENSG_A", "ENST_A", "rRNA");

        // overlap with A
        addGene(cache, "ENSG_B", "B", "5", 150, 250);
        addTranscript(cache, "ENSG_B", "ENST_B", "rRNA");

        // adjacent to merged [100..250] — touch should also merge
        addGene(cache, "ENSG_C", "C", "5", 251, 300);
        addTranscript(cache, "ENSG_C", "ENST_C", "rRNA");

        // disjoint
        addGene(cache, "ENSG_D", "D", "5", 1000, 1100);
        addTranscript(cache, "ENSG_D", "ENST_D", "rRNA");

        final RrnaRegionIndex index = RrnaRegionIndex.build(cache, RefGenomeVersion.V37);

        assertEquals(4, index.geneCount());
        // any point inside the merged 100..300 range hits
        assertTrue(index.overlaps("5", 200, 200));
        assertTrue(index.overlaps("5", 250, 251));
        assertTrue(index.overlaps("5", 300, 300));
        // gap between merged region and D
        assertFalse(index.overlaps("5", 301, 999));
        // D
        assertTrue(index.overlaps("5", 1050, 1050));
    }

    @Test
    public void testV38AddsAcrocentricPArms()
    {
        final EnsemblDataCache cache = new EnsemblDataCache("", RefGenomeVersion.V38);
        final RrnaRegionIndex index = RrnaRegionIndex.build(cache, RefGenomeVersion.V38);

        // empty cache + V38 -> biotype geneCount is 0 but acrocentric ranges are present
        assertEquals(0, index.geneCount());
        assertTrue(index.overlaps("chr13", 1, 1));
        assertTrue(index.overlaps("chr13", 15_999_999, 16_000_000));
        assertFalse(index.overlaps("chr13", 16_000_001, 16_000_001));
        assertTrue(index.overlaps("chr14", 5_000_000, 5_000_000));
        assertTrue(index.overlaps("chr15", 17_000_000, 17_000_000));
        assertTrue(index.overlaps("chr21", 11_999_999, 12_000_001));
        assertTrue(index.overlaps("chr22", 13_500_000, 13_500_000));
        // q-arms unaffected
        assertFalse(index.overlaps("chr13", 50_000_000, 50_000_001));
        // non-acrocentric chromosomes unaffected
        assertFalse(index.overlaps("chr1", 1_000_000, 1_000_000));
    }

    @Test
    public void testV38ExcludesChr22AndChrUnContigsByPrefix()
    {
        final EnsemblDataCache cache = new EnsemblDataCache("", RefGenomeVersion.V38);
        final RrnaRegionIndex index = RrnaRegionIndex.build(cache, RefGenomeVersion.V38);

        // wholesale exclusion regardless of position
        assertTrue(index.overlaps("chr22_KI270733v1_random", 1, 1));
        assertTrue(index.overlaps("chr22_KI270733v1_random", 100_000, 200_000));
        assertTrue(index.overlaps("chrUn_GL000220v1", 1, 1));
        assertTrue(index.overlaps("chrUn_GL000219v1", 5_000, 6_000));

        // main chr22 still falls under the acrocentric range, not the prefix rule
        assertTrue(index.overlaps("chr22", 1, 1));
        // chr22 q-arm not excluded
        assertFalse(index.overlaps("chr22", 30_000_000, 30_000_001));

        // chr13_/chr14_/chr15_/chr21_ random contigs are NOT excluded by prefix (only chr22_ / chrUn_)
        assertFalse(index.overlaps("chr21_KI270873v1_alt", 5_000, 6_000));
    }

    @Test
    public void testV37FallsBackToBiotypeOnly()
    {
        final EnsemblDataCache cache = new EnsemblDataCache("", RefGenomeVersion.V37);
        final RrnaRegionIndex index = RrnaRegionIndex.build(cache, RefGenomeVersion.V37);

        // V37 does NOT get acrocentric or prefix exclusions
        assertFalse(index.overlaps("chr13", 5_000_000, 5_000_000));
        assertFalse(index.overlaps("13", 5_000_000, 5_000_000));
        assertFalse(index.overlaps("chr22_KI270733v1_random", 1, 1));
        assertFalse(index.overlaps("chrUn_GL000220v1", 1, 1));
    }

    @Test
    public void testBuildSkipsGenesWithNoRrnaTranscript()
    {
        final EnsemblDataCache cache = new EnsemblDataCache("", RefGenomeVersion.V38);

        addGene(cache, "ENSG_MIXED", "MIXED", "1", 100, 200);
        // a gene with only non-rRNA transcripts is not added
        addTranscript(cache, "ENSG_MIXED", "ENST_A", "protein_coding");
        addTranscript(cache, "ENSG_MIXED", "ENST_B", "lncRNA");

        final RrnaRegionIndex index = RrnaRegionIndex.build(cache, RefGenomeVersion.V37);

        assertEquals(0, index.geneCount());
        assertFalse(index.overlaps("1", 150, 150));
    }

    private static void addGene(
            final EnsemblDataCache cache, final String geneId, final String geneName,
            final String chromosome, final int start, final int end)
    {
        final GeneData gene = new GeneData(geneId, geneName, chromosome, POS_STRAND, start, end, "");
        cache.getChrGeneDataMap().computeIfAbsent(chromosome, k -> new ArrayList<>()).add(gene);
    }

    private static void addTranscript(
            final EnsemblDataCache cache, final String geneId, final String transName, final String biotype)
    {
        final List<TranscriptData> transcripts = cache.getTranscriptDataMap()
                .computeIfAbsent(geneId, k -> new ArrayList<>());
        final TranscriptData transcript = new TranscriptData(
                transcripts.size() + 1, transName, geneId, false, POS_STRAND,
                0, 0, null, null, biotype, null);
        transcripts.add(transcript);
    }
}
