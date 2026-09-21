package com.hartwig.hmftools.viridian.integration.seq_align;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bwa.IBwaMemAligner;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;
import com.hartwig.hmftools.viridian.reference.ViralReference;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class ViralSequenceAlignerTest
{
    private static final OncologyGroup HPV_16 = new OncologyGroup("HPV 16");
    private static final OncologyGroup HPV_18 = new OncologyGroup("HPV 18");
    private static final ViralContig CONTIG_16 = new ViralContig("hpv16", 7906, "Human papillomavirus 16", HPV_16);
    private static final ViralContig CONTIG_18 = new ViralContig("hpv18", 7857, "Human papillomavirus 18", HPV_18);

    private static final String SEQUENCE = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    private static final String OTHER_SEQUENCE = "TTGCATTGCATTGCATTGCATTGCATTGCATT";

    @Test
    public void takesTheHighestScoringAlignment()
    {
        // The same sequence hits both contigs; only the better-scoring hit is kept.
        List<BwaMemAlignment> hits = List.of(
                alignment(0, 0, 999, 60, "32M", 1),
                alignment(0, 1, 499, 45, "32M", 4));

        assertEquals(
                new ViralSequenceAlignment(CONTIG_16, 1000, FORWARD, "32M", 60, 1, 32),
                align(SEQUENCE, hits));
    }

    // Near-identical contigs give a sequence equal scores, so the choice must not depend on the order BWA returns them.
    @Test
    public void breaksScoreTiesDeterministically()
    {
        List<BwaMemAlignment> ascending = List.of(
                alignment(0, 0, 999, 60, "32M", 1),
                alignment(0, 1, 499, 60, "32M", 1));
        List<BwaMemAlignment> descending = List.of(
                alignment(0, 1, 499, 60, "32M", 1),
                alignment(0, 0, 999, 60, "32M", 1));

        assertEquals(align(SEQUENCE, ascending), align(SEQUENCE, descending));
        assertEquals(CONTIG_16, align(SEQUENCE, descending).contig());
    }

    // A partly aligned sequence keeps the clip in its cigar, and its length stays that of the whole query, so how much
    // of the insert reached a virus can be recovered downstream.
    @Test
    public void reportsPartiallyAlignedSequence()
    {
        ViralSequenceAlignment alignment = align(SEQUENCE, List.of(alignment(0, 0, 199, 40, "12S20M", 0)));

        assertEquals("12S20M", alignment.cigar());
        assertEquals(32, alignment.sequenceLength());
    }

    @Test
    public void reportsReverseStrandAlignment()
    {
        assertEquals(
                new ViralSequenceAlignment(CONTIG_18, 300, REVERSE, "32M", 55, 2, 32),
                align(SEQUENCE, List.of(alignment(0x10, 1, 299, 55, "32M", 2))));
    }

    // A sequence matching no contig yields no alignment, and the results stay aligned with the input order.
    @Test
    public void returnsNullWhereNothingAligned()
    {
        ViralSequenceAligner aligner = aligner(Map.of(
                OTHER_SEQUENCE, List.of(noHit()),
                SEQUENCE, List.of(alignment(0, 0, 999, 60, "32M", 1))));

        List<ViralSequenceAlignment> alignments = aligner.alignAll(List.of(OTHER_SEQUENCE, SEQUENCE));

        assertNull(alignments.get(0));
        assertEquals(CONTIG_16, alignments.get(1).contig());
    }

    private static ViralSequenceAlignment align(String sequence, List<BwaMemAlignment> hits)
    {
        return aligner(Map.of(sequence, hits)).alignAll(List.of(sequence)).get(0);
    }

    private static ViralSequenceAligner aligner(Map<String, List<BwaMemAlignment>> hits)
    {
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord(CONTIG_16.name(), CONTIG_16.length()),
                new SAMSequenceRecord(CONTIG_18.name(), CONTIG_18.length())));
        ViralReference reference = new ViralReference(List.of(CONTIG_16, CONTIG_18), dictionary);

        return new ViralSequenceAligner(new FakeAligner(hits), reference);
    }

    private static BwaMemAlignment alignment(int samFlag, int refId, int refStart, int score, String cigar, int editDistance)
    {
        return new BwaMemAlignment(
                samFlag, refId, refStart, refStart + 32, 0, 32, 60, editDistance, score, 0, cigar, null, null, -1, -1, 0);
    }

    private static BwaMemAlignment noHit()
    {
        return new BwaMemAlignment(0x4, -1, -1, -1, 0, 0, 0, 0, 0, 0, "*", null, null, -1, -1, 0);
    }

    private record FakeAligner(Map<String, List<BwaMemAlignment>> hitsBySequence) implements IBwaMemAligner
    {
        @Override
        public List<List<BwaMemAlignment>> alignSequences(List<byte[]> sequences)
        {
            return sequences.stream().map(sequence -> hitsBySequence.get(new String(sequence))).toList();
        }
    }
}
