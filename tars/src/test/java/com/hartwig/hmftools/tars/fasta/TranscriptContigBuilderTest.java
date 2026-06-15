package com.hartwig.hmftools.tars.fasta;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class TranscriptContigBuilderTest
{
    private static final String GENE_ID = "ENSG00000171862";
    private static final String GENE_NAME = "PTEN";
    private static final String TRANS_NAME = "ENST00000399770";

    @Test
    public void testThreeExonForwardStrandTranscript()
    {
        final TranscriptData transcript = transcript();
        addExon(transcript, 100, 199, 1);
        addExon(transcript, 300, 399, 2);
        addExon(transcript, 500, 549, 3);

        final TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        final TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene(POS_STRAND), transcript);

        assertNotNull(result);
        assertEquals(GENE_ID, result.geneId());
        assertEquals(GENE_NAME, result.geneName());
        assertEquals(TRANS_NAME, result.transName());
        assertEquals(CHR_1, result.chromosome());

        assertEquals(3, result.exonSpans().size());
        assertEquals(new BaseRegion(100, 199), result.exonSpans().get(0));
        assertEquals(new BaseRegion(300, 399), result.exonSpans().get(1));
        assertEquals(new BaseRegion(500, 549), result.exonSpans().get(2));

        assertEquals(100 + 100 + 50, result.sequence().length());
    }

    @Test
    public void testNegativeStrandTranscriptStoresExonsInGenomicOrder()
    {
        // Rank-order for a negative-strand transcript is genomic-descending; builder must re-sort to genomic-ascending.
        final TranscriptData transcript = transcript();
        addExon(transcript, 500, 549, 1); // 5' end of mRNA, but genomic-high
        addExon(transcript, 300, 399, 2);
        addExon(transcript, 100, 199, 3);

        final TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        final TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene(NEG_STRAND), transcript);

        assertNotNull(result);
        assertEquals(100, result.exonSpans().get(0).start());
        assertEquals(300, result.exonSpans().get(1).start());
        assertEquals(500, result.exonSpans().get(2).start());
    }

    @Test
    public void testEmittedSequenceMatchesExonConcatenation()
    {
        final TranscriptData transcript = transcript();
        addExon(transcript, 100, 199, 1);
        addExon(transcript, 300, 399, 2);
        addExon(transcript, 500, 549, 3);

        final MockRefGenome ref = buildMockRef();
        final TranscriptContigBuilder builder = new TranscriptContigBuilder(ref);
        final TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene(POS_STRAND), transcript);

        assertNotNull(result);

        final String expected = ref.getBaseString(CHR_1, 100, 199)
                + ref.getBaseString(CHR_1, 300, 399)
                + ref.getBaseString(CHR_1, 500, 549);
        assertEquals(expected, result.sequence());
    }

    @Test
    public void testMultipleTranscriptsEachGetTheirOwnContig()
    {
        final TranscriptData transA = transcript(1, "ENST00000001", true);
        addExon(transA, 100, 199, 1);
        addExon(transA, 300, 399, 2);

        final TranscriptData transB = transcript(2, "ENST00000002", false);
        addExon(transB, 100, 199, 1);
        addExon(transB, 500, 549, 2);

        final TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        final TranscriptContigBuilder.TranscriptContigResult resultA = builder.build(gene(POS_STRAND), transA);
        final TranscriptContigBuilder.TranscriptContigResult resultB = builder.build(gene(POS_STRAND), transB);

        assertNotNull(resultA);
        assertNotNull(resultB);

        assertEquals("ENST00000001", resultA.transName());
        assertEquals("ENST00000002", resultB.transName());

        assertEquals(2, resultA.exonSpans().size());
        assertEquals(2, resultB.exonSpans().size());

        assertEquals(200, resultA.sequence().length());
        assertEquals(150, resultB.sequence().length());
    }

    @Test
    public void testTranscriptWithNoExonsReturnsNull()
    {
        final TranscriptData empty = transcript();
        final TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        assertNull(builder.build(gene(POS_STRAND), empty));
    }

    @Test
    public void testSingleExonTranscriptHasNoIntrons()
    {
        final TranscriptData transcript = transcript();
        addExon(transcript, 100, 200, 1);

        final TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        final TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene(POS_STRAND), transcript);

        assertNotNull(result);
        assertEquals(1, result.exonSpans().size());
        assertEquals(101, result.sequence().length());
    }

    private static GeneData gene(final byte strand)
    {
        return createEnsemblGeneData(GENE_ID, GENE_NAME, CHR_1, strand, 100, 549);
    }

    // Local builders are retained because GeneTestUtils.createTransExons only models uniform exon length and
    // genomic-ascending rank, whereas these tests need non-uniform exons (100/100/50bp) and a negative-strand
    // transcript whose ranks descend genomically.
    private static TranscriptData transcript()
    {
        return transcript(1, TRANS_NAME, true);
    }

    private static TranscriptData transcript(final int transId, final String transName, final boolean canonical)
    {
        return new TranscriptData(transId, transName, GENE_ID, canonical, POS_STRAND, 0, 0, null, null, "protein_coding", null);
    }

    private static void addExon(final TranscriptData transcript, final int start, final int end, final int rank)
    {
        transcript.exons().add(new ExonData(transcript.TransId, start, end, rank, -1, -1));
    }

    private static MockRefGenome buildMockRef()
    {
        // position-distinguishable bases so sequence-content assertions actually catch bad concatenations
        MockRefGenome refGenome = new MockRefGenome(true);
        char[] alphabet = { 'A', 'C', 'G', 'T' };
        StringBuilder bases = new StringBuilder();
        bases.append('N'); // pad index 0
        for(int i = 1; i <= 1000; ++i)
            bases.append(alphabet[i % 4]);
        refGenome.RefGenomeMap.put(CHR_1, bases.toString());
        return refGenome;
    }
}
