package com.hartwig.hmftools.geneutils.probequality;

import static java.lang.Math.min;

import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.util.Comparator;
import java.util.List;
import java.util.function.Supplier;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

// Evaluates the off-target risk of a probe given its alignments from BWA.
public class ProbeQualityModel
{
    private final BwaMemAligner mAligner;
    // Desired length of the probes.
    // All probes tested with this model must have this length (BWA-MEM config relies on it).
    private final int mTargetProbeLength;
    // Alignment score must exceed this to count towards the risk score.
    private final int mMatchScoreThreshold;
    // Amount that 1 alignment match counts towards the risk score.
    // E.g. value of 10 means an alignment with score = mMatchScoreThreshold contributes 10 risk score points.
    private final int mMatchScoreOffset;

    public ProbeQualityModel(final Supplier<BwaMemAligner> alignerFactory, final int mTargetProbeLength, final int mMatchScoreThreshold,
            final int mMatchScoreOffset)
    {
        mAligner = alignerFactory.get();
        if(mTargetProbeLength < 1)
        {
            throw new IllegalArgumentException("targetProbeLength must be >= 1");
        }
        this.mTargetProbeLength = mTargetProbeLength;
        this.mMatchScoreThreshold = mMatchScoreThreshold;
        if(mMatchScoreOffset < 0)
        {
            // < 0 will break the maths.
            throw new IllegalArgumentException("matchScoreOffset must be >= 0");
        }
        this.mMatchScoreOffset = mMatchScoreOffset;

        setBwaMemParams(mAligner);
        logBwaMemParams(mAligner);
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

    // Compute probe qualities for a list of probes.
    public List<Result> compute(final List<byte[]> probes)
    {
        probes.forEach(probe ->
        {
            if(probe.length != mTargetProbeLength)
            {
                throw new IllegalArgumentException("Probe sequence length doesn't match configured probe length");
            }
        });

        GU_LOGGER.debug("Running BWA-MEM alignment");
        List<List<BwaMemAlignment>> alignments = mAligner.alignSeqs(probes);
        if(alignments.size() != probes.size())
        {
            // Presumably this shouldn't occur, but we'll check to give a nicer error just in case.
            throw new RuntimeException("Alignment failed");
        }
        GU_LOGGER.debug("Running risk model");
        return alignments.stream().map(this::computeFromAlignments).toList();
    }

    public Result computeFromAlignments(final List<BwaMemAlignment> alignments)
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

    private void setBwaMemParams(BwaMemAligner aligner)
    {
        // Ensure we can find alignments fitting our parameters.
        aligner.setMinSeedLengthOption(min(min(19, mMatchScoreThreshold + 10), mTargetProbeLength / 2));
        // Output many alignments per query.
        aligner.setFlagOption(aligner.getFlagOption() | BwaMemAligner.MEM_F_ALL);
        aligner.setOutputScoreThresholdOption(mMatchScoreThreshold);
        // Don't prune seeds with many occurrences in the genome. This is a key performance tuning parameter.
        aligner.setMaxMemIntvOption(2000);
        aligner.setMaxSeedOccurencesOption(10000);
        // Other minor params to encourage more alignments to be found.
        aligner.setBandwidthOption(mTargetProbeLength);
        aligner.setSplitFactorOption(0.5f);
        aligner.setDropRatioOption(0.1f);
    }

    private static void logBwaMemParams(BwaMemAligner aligner)
    {
        GU_LOGGER.debug("BWA-MEM options:");
        GU_LOGGER.debug("  MinSeedLength: {}", aligner.getMinSeedLengthOption());
        GU_LOGGER.debug("  SplitFactor: {}", aligner.getSplitFactorOption());
        GU_LOGGER.debug("  SplitWidth: {}", aligner.getSplitWidthOption());
        GU_LOGGER.debug("  MaxSeedOccurrences: {}", aligner.getMaxSeedOccurencesOption());
        GU_LOGGER.debug("  MaxMemOccurrences: {}", aligner.getMaxMemIntvOption());
        GU_LOGGER.debug("  DropRatio: {}", aligner.getDropRatioOption());
        GU_LOGGER.debug("  Match: {}", aligner.getMatchScoreOption());
        GU_LOGGER.debug("  Mismatch: {}", aligner.getMismatchPenaltyOption());
        GU_LOGGER.debug("  IGapOpen: {}", aligner.getIGapOpenPenaltyOption());
        GU_LOGGER.debug("  IGapExtend: {}", aligner.getIGapExtendPenaltyOption());
        GU_LOGGER.debug("  DGapOpen: {}", aligner.getDGapOpenPenaltyOption());
        GU_LOGGER.debug("  DGapExtend: {}", aligner.getDGapExtendPenaltyOption());
        GU_LOGGER.debug("  Clip3: {}", aligner.getClip3PenaltyOption());
        GU_LOGGER.debug("  Clip5: {}", aligner.getClip5PenaltyOption());
        GU_LOGGER.debug("  Bandwidth: {}", aligner.getBandwidthOption());
        GU_LOGGER.debug("  ZDrop: {}", aligner.getZDropOption());
        GU_LOGGER.debug("  OutputScoreThreshold: {}", aligner.getOutputScoreThresholdOption());
        GU_LOGGER.debug("  Flags: {}", aligner.getFlagOption());
        GU_LOGGER.debug("  Threads: {}", aligner.getNThreadsOption());
    }
}
