package com.hartwig.hmftools.panelbuilder.probequality;

import static java.lang.Math.min;

import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.bwa.BwaMemAlignParams;
import com.hartwig.hmftools.common.bwa.BwaMemAlignerConfig;
import com.hartwig.hmftools.common.bwa.IBwaMemAligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

// Evaluates the off-target risk of a probe given its alignments from BWA.
public class ProbeQualityModel
{
    private final IBwaMemAligner mAligner;
    // Desired length of the probes.
    // All probes tested with this model must have this length (BWA-MEM config relies on it).
    private final int mTargetProbeLength;
    // Alignment score must exceed this to count towards the risk score.
    private final int mMatchScoreThreshold;
    // Amount that one alignment match counts towards the risk score.
    // E.g. value of 10 means an alignment with score = mMatchScoreThreshold contributes 10 risk score points.
    private final int mMatchScoreOffset;

    // Use in production; the constructor is for injecting a test aligner.
    public static ProbeQualityModel create(
            final String bwaIndexImage, final int threads, final int targetProbeLength, final int matchScoreThreshold,
            final int matchScoreOffset)
    {
        BwaMemAlignerConfig config = new BwaMemAlignerConfig(
                bwaIndexImage, alignParams(targetProbeLength, matchScoreThreshold), true, threads, null);
        return new ProbeQualityModel(new BwaMemAligner(config), targetProbeLength, matchScoreThreshold, matchScoreOffset);
    }

    public ProbeQualityModel(final IBwaMemAligner aligner, final int targetProbeLength, final int matchScoreThreshold,
            final int matchScoreOffset)
    {
        mAligner = aligner;
        if(targetProbeLength < 1)
        {
            throw new IllegalArgumentException("targetProbeLength must be >= 1");
        }
        mTargetProbeLength = targetProbeLength;
        mMatchScoreThreshold = matchScoreThreshold;
        if(matchScoreOffset < 0)
        {
            // < 0 will break the maths.
            throw new IllegalArgumentException("matchScoreOffset must be >= 0");
        }
        mMatchScoreOffset = matchScoreOffset;
    }

    // BWA-MEM alignment parameters tuned to find the many low-scoring off-target matches this model relies on.
    private static BwaMemAlignParams alignParams(int targetProbeLength, int matchScoreThreshold)
    {
        return BwaMemAlignParams.DEFAULT
                .withSeedLengthMin(min(19, targetProbeLength / 2))
                // Don't prune seeds with many genome occurrences (key performance parameter).
                .withSeed3MaxOccurrence(2000)
                .withMemMaxOccurrence(2000)
                // Other minor params to encourage more alignments to be found.
                .withMemReseedFactor(0.5f)
                .withChainOverlapFactor(0.1f)
                .withBandWidth(targetProbeLength)
                .withMinAlignScore(matchScoreThreshold);
    }

    public record Result(
            // Range: [0, 1].
            // Roughly reciprocal to the number of exact match alignments. E.g. 1 means no off-target, 0.5 means 1 off-target.
            double qualityScore,
            // Range: [0, +inf].
            // Higher is higher chance of off-target match.
            long riskScore,
            // Number of alignments contributing to the risk score (i.e. above the threshold).
            int offTargetCount,
            // Raw sum of alignment scores of alignments contributing to the risk score.
            long offTargetScoreSum
    )
    {
    }

    public List<Result> computeFromSeqString(final List<String> probes)
    {
        List<byte[]> probeBytes = probes.stream().map(String::getBytes).toList();
        return computeFromSeqBytes(probeBytes);
    }

    // Compute probe qualities for a list of probes.
    public List<Result> computeFromSeqBytes(final List<byte[]> probes)
    {
        probes.forEach(probe ->
        {
            if(probe.length != mTargetProbeLength)
            {
                throw new IllegalArgumentException("Probe sequence length doesn't match configured probe length");
            }
        });

        List<List<BwaMemAlignment>> alignments = mAligner.alignSequences(probes);
        return alignments.stream().map(this::computeFromAlignments).toList();
    }

    private Result computeFromAlignments(final List<BwaMemAlignment> alignments)
    {
        // Order by best match first.
        alignments.sort(Comparator.comparing(BwaMemAlignment::getAlignerScore, Comparator.reverseOrder()));
        // First alignment which is assumed to be the on-target exact match.
        int targetScore = alignments.get(0).getAlignerScore();
        List<Integer> offTarget = alignments.stream()
                // Drop the first alignment which is assumed to be the on-target exact match.
                .skip(1)
                // Only need the alignment scores. The alignment score from BWA-MEM is effectively a similarity score.
                .map(BwaMemAlignment::getAlignerScore)
                // Keep only alignments with score above the configured threshold.
                .takeWhile(score -> score >= mMatchScoreThreshold)
                .toList();
        int offTargetCount = offTarget.size();
        long offTargetScoreSum = offTarget.stream().mapToLong(s -> s).sum();
        long riskScore = offTargetScoreSum - (long) offTargetCount * (mMatchScoreThreshold - mMatchScoreOffset);
        // This value is the approximate equivalent total count of exact match alignments.
        double effectiveOffTargetMatchLength = (double) riskScore / (targetScore - mMatchScoreThreshold + mMatchScoreOffset);
        double qualityScore = 1 / (1 + effectiveOffTargetMatchLength);
        return new Result(qualityScore, riskScore, offTargetCount, offTargetScoreSum);
    }
}
