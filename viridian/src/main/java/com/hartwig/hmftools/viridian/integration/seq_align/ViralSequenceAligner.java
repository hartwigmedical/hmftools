package com.hartwig.hmftools.viridian.integration.seq_align;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATION_ALIGN_SCORE_MIN;

import java.util.Comparator;
import java.util.List;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.bwa.BwaMemAlignParams;
import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.bwa.BwaMemAlignerConfig;
import com.hartwig.hmftools.common.bwa.IBwaMemAligner;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.viridian.reference.ViralContig;
import com.hartwig.hmftools.viridian.reference.ViralReference;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFlag;

// Aligns candidate host integration variant sequences to the viral reference.
// Produces the only best alignment across the entire viral reference - for data gathering purposes.
public class ViralSequenceAligner
{
    private final IBwaMemAligner mAligner;
    private final ViralReference mReference;

    public static ViralSequenceAligner create(ViralReference reference, String bwaIndexImage, int threads)
    {
        BwaMemAlignParams params = BwaMemAlignParams.DEFAULT.withMinAlignScore(INTEGRATION_ALIGN_SCORE_MIN);
        BwaMemAlignerConfig alignerConfig = new BwaMemAlignerConfig(bwaIndexImage, params, false, threads, null);
        BwaMemAligner aligner = new BwaMemAligner(alignerConfig);
        return new ViralSequenceAligner(aligner, reference);
    }

    ViralSequenceAligner(IBwaMemAligner aligner, ViralReference reference)
    {
        mAligner = aligner;
        mReference = reference;
    }

    // One result per input sequence, in the same order. An entry is null where the sequence is aligned to no contig.
    public List<ViralSequenceAlignment> alignAll(List<String> sequences)
    {
        List<byte[]> queries = sequences.stream().map(String::getBytes).toList();
        List<List<BwaMemAlignment>> alignments = mAligner.alignSequences(queries);

        return IntStream.range(0, sequences.size())
                .mapToObj(i -> selectBestAlignment(alignments.get(i), sequences.get(i).length()))
                .toList();
    }

    @Nullable
    private ViralSequenceAlignment selectBestAlignment(List<BwaMemAlignment> alignments, int sequenceLength)
    {
        return alignments.stream()
                .filter(alignment -> alignment.getRefId() >= 0)
                .min(BEST_ALIGNMENT_FIRST)
                .map(alignment -> toViralSequenceAlignment(alignment, sequenceLength))
                .orElse(null);
    }

    private ViralSequenceAlignment toViralSequenceAlignment(BwaMemAlignment alignment, int sequenceLength)
    {
        String contigName = mReference.sequenceDictionary().getSequence(alignment.getRefId()).getSequenceName();
        ViralContig viralContig = mReference.contig(contigName);
        // BWA reference positions are 0-based, so convert to 1-based.
        int position = alignment.getRefStart() + 1;
        Orientation orientation = SamRecordUtils.isFlagSet(alignment.getSamFlag(), SAMFlag.READ_REVERSE_STRAND)
                ? Orientation.REVERSE : Orientation.FORWARD;
        return new ViralSequenceAlignment(
                viralContig,
                position,
                orientation,
                cigarFromStr(alignment.getCigar()),
                alignment.getAlignerScore(),
                alignment.getNMismatches(),
                sequenceLength);
    }

    private static final Comparator<BwaMemAlignment> BEST_ALIGNMENT_FIRST = Comparator
            .<BwaMemAlignment>comparingInt(a -> -a.getAlignerScore())
            // Deterministic tiebreakers.
            .thenComparing(BwaMemAlignment::getRefId)
            .thenComparing(BwaMemAlignment::getRefStart);
}
