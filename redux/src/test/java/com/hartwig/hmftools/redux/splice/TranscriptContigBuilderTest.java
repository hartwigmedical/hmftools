package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

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
        // exons stored in transcribed (Rank) order, which equals genomic-ascending for a forward-strand gene
        TranscriptData transcript = transcript();
        addExon(transcript, 100, 199, 1);
        addExon(transcript, 300, 399, 2);
        addExon(transcript, 500, 549, 3);

        TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene(POS_STRAND), transcript);

        assertNotNull(result);
        assertEquals("ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME, result.Entry.contigName());
        assertEquals(GENE_ID, result.Entry.geneId());
        assertEquals(GENE_NAME, result.Entry.geneName());
        assertEquals(TRANS_NAME, result.Entry.transName());

        assertEquals(3, result.Entry.exonSpans().size());
        assertEquals(new BaseRegion(100, 199), result.Entry.exonSpans().get(0));
        assertEquals(new BaseRegion(300, 399), result.Entry.exonSpans().get(1));
        assertEquals(new BaseRegion(500, 549), result.Entry.exonSpans().get(2));

        // sequence length equals total exonic bases (no introns)
        assertEquals(100 + 100 + 50, result.Sequence.length());
        assertEquals(result.Entry.contigLength(), result.Sequence.length());
    }

    @Test
    public void testNegativeStrandTranscriptStoresExonsInGenomicOrder()
    {
        // Rank-ordered exons for a negative-strand transcript run genomic-descending. The builder must re-sort to
        // genomic-ascending so the contig sequence mirrors the forward strand.
        TranscriptData transcript = transcript();
        addExon(transcript, 500, 549, 1); // 5' end of mRNA, but genomic-high
        addExon(transcript, 300, 399, 2);
        addExon(transcript, 100, 199, 3);

        TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene(NEG_STRAND), transcript);

        assertNotNull(result);
        assertEquals(100, result.Entry.exonSpans().get(0).start());
        assertEquals(300, result.Entry.exonSpans().get(1).start());
        assertEquals(500, result.Entry.exonSpans().get(2).start());
    }

    @Test
    public void testEmittedSequenceMatchesExonConcatenation()
    {
        TranscriptData transcript = transcript();
        addExon(transcript, 100, 199, 1);
        addExon(transcript, 300, 399, 2);
        addExon(transcript, 500, 549, 3);

        MockRefGenome ref = buildMockRef();
        TranscriptContigBuilder builder = new TranscriptContigBuilder(ref);
        TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene(POS_STRAND), transcript);

        assertNotNull(result);

        String expected = ref.getBaseString(CHR_1, 100, 199)
                + ref.getBaseString(CHR_1, 300, 399)
                + ref.getBaseString(CHR_1, 500, 549);
        assertEquals(expected, result.Sequence);
    }

    @Test
    public void testMultipleTranscriptsEachGetTheirOwnContig()
    {
        // simulate two isoforms of the same gene with different exon sets — each should be built independently
        // (caller iterates transcripts; the builder is per-transcript)
        TranscriptData transA = new TranscriptData(1, "ENST00000001", GENE_ID, true, POS_STRAND, 0, 0, null, null, "protein_coding", null);
        addExon(transA, 100, 199, 1);
        addExon(transA, 300, 399, 2);

        TranscriptData transB = new TranscriptData(2, "ENST00000002", GENE_ID, false, POS_STRAND, 0, 0, null, null, "protein_coding", null);
        addExon(transB, 100, 199, 1);
        addExon(transB, 500, 549, 2);

        TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        TranscriptContigBuilder.TranscriptContigResult resultA = builder.build(gene(POS_STRAND), transA);
        TranscriptContigBuilder.TranscriptContigResult resultB = builder.build(gene(POS_STRAND), transB);

        assertNotNull(resultA);
        assertNotNull(resultB);

        assertEquals("ens" + GENE_ID + "_" + GENE_NAME + "_ENST00000001", resultA.Entry.contigName());
        assertEquals("ens" + GENE_ID + "_" + GENE_NAME + "_ENST00000002", resultB.Entry.contigName());

        assertEquals(2, resultA.Entry.exonSpans().size());
        assertEquals(2, resultB.Entry.exonSpans().size());

        // contig A spans 100-199 + 300-399, contig B spans 100-199 + 500-549. Different lengths confirm they're independent.
        assertEquals(200, resultA.Entry.contigLength());
        assertEquals(150, resultB.Entry.contigLength());
    }

    @Test
    public void testTranscriptWithNoExonsReturnsNull()
    {
        TranscriptData empty = transcript();
        TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        assertNull(builder.build(gene(POS_STRAND), empty));
    }

    @Test
    public void testSingleExonTranscriptHasNoIntrons()
    {
        TranscriptData transcript = transcript();
        addExon(transcript, 100, 200, 1);

        TranscriptContigBuilder builder = new TranscriptContigBuilder(buildMockRef());
        TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene(POS_STRAND), transcript);

        assertNotNull(result);
        assertEquals(1, result.Entry.exonSpans().size());
        assertEquals(101, result.Sequence.length());
    }

    private static GeneData gene(final byte strand)
    {
        return new GeneData(GENE_ID, GENE_NAME, CHR_1, strand, 100, 549, "");
    }

    private static TranscriptData transcript()
    {
        return new TranscriptData(1, TRANS_NAME, GENE_ID, true, POS_STRAND, 0, 0, null, null, "protein_coding", null);
    }

    private static void addExon(final TranscriptData transcript, int start, int end, int rank)
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
