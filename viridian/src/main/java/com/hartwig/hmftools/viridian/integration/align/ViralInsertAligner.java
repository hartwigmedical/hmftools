package com.hartwig.hmftools.viridian.integration.align;

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
import com.hartwig.hmftools.viridian.reference.VirusReference;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFlag;

// Aligns candidate integration SV insert sequences to the virus reference.
// Produces the only best alignment across the entire virus reference.
// TODO: do we need to constrain the alignment to detected virus genomes only?
public class ViralInsertAligner
{
    private final IBwaMemAligner mAligner;
    private final VirusReference mReference;

    public static ViralInsertAligner create(VirusReference reference, String bwaIndexImage, int threads)
    {
        BwaMemAlignParams params = BwaMemAlignParams.DEFAULT.withMinAlignScore(INTEGRATION_ALIGN_SCORE_MIN);
        BwaMemAlignerConfig alignerConfig = new BwaMemAlignerConfig(bwaIndexImage, params, false, threads, null);
        BwaMemAligner aligner = new BwaMemAligner(alignerConfig);
        return new ViralInsertAligner(aligner, reference);
    }

    ViralInsertAligner(IBwaMemAligner aligner, VirusReference reference)
    {
        mAligner = aligner;
        mReference = reference;
    }

    // One result per input sequence, in the same order. An entry is null where the sequence is aligned to no contig.
    public List<ViralInsertAlignment> align(List<String> sequences)
    {
        List<byte[]> queries = sequences.stream().map(String::getBytes).toList();
        List<List<BwaMemAlignment>> alignments = mAligner.alignSequences(queries);

        return IntStream.range(0, sequences.size())
                .mapToObj(i -> selectBestAlignment(alignments.get(i), sequences.get(i).length()))
                .toList();
    }

    @Nullable
    private ViralInsertAlignment selectBestAlignment(List<BwaMemAlignment> alignments, int sequenceLength)
    {
        return alignments.stream()
                .filter(alignment -> alignment.getRefId() >= 0)
                .min(BEST_ALIGNMENT_FIRST)
                .map(alignment -> toViralInsertAlignment(alignment, sequenceLength))
                .orElse(null);
    }

    private ViralInsertAlignment toViralInsertAlignment(BwaMemAlignment alignment, int sequenceLength)
    {
        String contigName = mReference.sequenceDictionary().getSequence(alignment.getRefId()).getSequenceName();
        ViralContig viralContig = mReference.contig(contigName);
        // BWA reference positions are 0-based, so convert to 1-based.
        int position = alignment.getRefStart() + 1;
        Orientation orientation = SamRecordUtils.isFlagSet(alignment.getSamFlag(), SAMFlag.READ_REVERSE_STRAND)
                ? Orientation.REVERSE : Orientation.FORWARD;
        return new ViralInsertAlignment(
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
