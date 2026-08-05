package com.hartwig.hmftools.panelbuilder.probequality;

import static java.lang.Math.min;

import static com.hartwig.hmftools.panelbuilder.probequality.Utils.partitionStream;

import java.util.List;
import java.util.Map;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
    // Amount that one alignment match counts towards the risk score.
    // E.g. value of 10 means an alignment with score = mMatchScoreThreshold contributes 10 risk score points.
    private final int mMatchScoreOffset;
    // BWA-MEM per-base match reward; a full-length perfect alignment scores mTargetProbeLength * mMatchScore.
    private final int mMatchScore;
    // Maps a BWA alignment's reference contig index (getRefId) to its chromosome name, so an alignment can be tested against a probe's source
    // regions. Built from the reference sequence dictionary, whose order matches the BWA index.
    private final Map<Integer, String> mRefIdToChromosome;

    private static final Logger LOGGER = LogManager.getLogger(ProbeQualityModel.class);

    public ProbeQualityModel(final Supplier<BwaMemAligner> alignerFactory, final int mTargetProbeLength, final int mMatchScoreThreshold,
            final int mMatchScoreOffset, final Map<Integer, String> refIdToChromosome)
    {
        mAligner = alignerFactory.get();
        mRefIdToChromosome = refIdToChromosome;
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
        mMatchScore = mAligner.getMatchScoreOption();
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

    // Builds the refId -> chromosome map the model needs, from a reference genome sequence dictionary. The dictionary order matches the BWA
    // index contig order, so an alignment's getRefId() indexes into it.
    public static Map<Integer, String> buildRefIdToChromosome(final SAMSequenceDictionary refGenomeDictionary)
    {
        return refGenomeDictionary.getSequences().stream()
                .collect(Collectors.toMap(SAMSequenceRecord::getSequenceIndex, SAMSequenceRecord::getSequenceName));
    }

    public List<Result> computeFromSeqString(final List<String> probes, final List<List<ChrBaseRegion>> probeSourceRegions)
    {
        List<byte[]> probeBytes = probes.stream().map(String::getBytes).toList();
        return computeFromSeqBytes(probeBytes, probeSourceRegions);
    }

    // Compute probe qualities for a list of probes. probeSourceRegions[i] is the reference region(s) probe i is built from (its intended
    // capture loci): one region for a normal reference-contiguous probe, several for a constructed probe (each fragment of a variant/SV probe
    // or RNA spliced probe). Alignments landing on those regions are the on-target captures and are excluded from the off-target risk; the
    // rest count as off-target. The on-target quality is normalised against a theoretical full-length match.
    public List<Result> computeFromSeqBytes(final List<byte[]> probes, final List<List<ChrBaseRegion>> probeSourceRegions)
    {
        if(probes.size() != probeSourceRegions.size())
        {
            throw new IllegalArgumentException("probes and probeSourceRegions must be the same size");
        }
        probes.forEach(probe ->
        {
            if(probe.length != mTargetProbeLength)
            {
                throw new IllegalArgumentException("Probe sequence length doesn't match configured probe length");
            }
        });

        List<List<BwaMemAlignment>> alignments = runAlignment(probes, 0);
        return IntStream.range(0, alignments.size())
                .mapToObj(i -> computeFromAlignments(alignments.get(i), probeSourceRegions.get(i)))
                .toList();
    }

    private List<List<BwaMemAlignment>> runAlignment(List<byte[]> probes, int batchSize)
    {
        if(batchSize <= 0)
        {
            batchSize = probes.size();
        }
        try
        {
            Stream<List<byte[]>> batches = batchSize < probes.size() ? partitionStream(probes.stream(), batchSize) : Stream.of(probes);
            return batches
                    .flatMap(batch ->
                    {
                        LOGGER.trace("Running BWA-MEM alignment on {} queries", batch.size());
                        List<List<BwaMemAlignment>> batchAlignments = mAligner.alignSeqs(batch);
                        if(batchAlignments.size() != batch.size())
                        {
                            // Presumably, this shouldn't occur, but we'll check to give a nicer error just in case.
                            throw new RuntimeException("Alignment failed");
                        }
                        return batchAlignments.stream();
                    })
                    .toList();
        }
        catch(IllegalArgumentException e)
        {
            // Silly issue with the BWA-MEM JNI wrapper where it may try to allocate a ByteBuffer larger than is allowed (which is
            // subsequently a silly issue with ByteBuffer being limited to 2^31 bytes).
            // We don't have much control over this issue since the number of alignments produced can't be predicted or tightly controlled.
            // Further, modifying the wrapper code is a pain.
            // So instead we will catch the error and shrink the sequence set, hopefully reducing the allocation to an acceptable size.
            int newBatchSize = batchSize / 10;
            if(newBatchSize >= 1)
            {
                LOGGER.warn(
                        "Aligning sequences with batch size {} failed, trying again with batch size {}",
                        batchSize, newBatchSize);
                return runAlignment(probes, newBatchSize);
            }
            else
            {
                throw e;
            }
        }
    }

    private Result computeFromAlignments(final List<BwaMemAlignment> alignments, final List<ChrBaseRegion> sourceRegions)
    {
        // Normalise against a theoretical perfect full-length match. The probe's intended target is its source region(s): for a normal probe,
        // the single region it was taken from; for a constructed probe (variant/SV, RNA spliced), each fragment's source locus - all of which
        // it is designed to capture. Alignments landing on those regions are the on-target captures and are excluded; every other alignment
        // above the threshold is off-target. (Previously the single top alignment was assumed to be a full-length on-target match and dropped,
        // and its score used as the normaliser - both wrong for constructed sequences, whose best hit is only a partial fragment.)
        int targetScore = mTargetProbeLength * mMatchScore;
        List<Integer> offTarget = alignments.stream()
                .filter(alignment -> !isOnTarget(alignment, sourceRegions))
                // Only need the alignment scores. The alignment score from BWA-MEM is effectively a similarity score.
                .map(BwaMemAlignment::getAlignerScore)
                .filter(score -> score >= mMatchScoreThreshold)
                .toList();
        int offTargetCount = offTarget.size();
        long offTargetScoreSum = offTarget.stream().mapToLong(s -> s).sum();
        long riskScore = offTargetScoreSum - (long) offTargetCount * (mMatchScoreThreshold - mMatchScoreOffset);
        // This value is the approximate equivalent total count of exact match alignments.
        double effectiveOffTargetMatchLength = (double) riskScore / (targetScore - mMatchScoreThreshold + mMatchScoreOffset);
        double qualityScore = 1 / (1 + effectiveOffTargetMatchLength);
        return new Result(qualityScore, riskScore, offTargetCount, offTargetScoreSum);
    }

    // Whether an alignment lands on one of the probe's source regions (an intended on-target capture). The alignment's reference contig index
    // is resolved to a chromosome name via the reference sequence dictionary, then compared to the source regions.
    private boolean isOnTarget(final BwaMemAlignment alignment, final List<ChrBaseRegion> sourceRegions)
    {
        String chromosome = mRefIdToChromosome.get(alignment.getRefId());
        if(chromosome == null)
        {
            return false;
        }
        int alignmentStart = alignment.getRefStart() + 1;   // BWA reference positions are 0-based; source regions are 1-based.
        int alignmentEnd = alignment.getRefEnd();
        for(ChrBaseRegion region : sourceRegions)
        {
            if(chromosome.equals(region.chromosome()) && alignmentStart <= region.end() && alignmentEnd >= region.start())
            {
                return true;
            }
        }
        return false;
    }

    private void setBwaMemParams(BwaMemAligner aligner)
    {
        // Ensure we can find alignments fitting our parameters.
        aligner.setMinSeedLengthOption(min(19, mTargetProbeLength / 2));
        // Output many alignments per query.
        aligner.setFlagOption(aligner.getFlagOption() | BwaMemAligner.MEM_F_ALL);
        aligner.setOutputScoreThresholdOption(mMatchScoreThreshold);
        // Don't prune seeds with many occurrences in the genome. This is a key performance tuning parameter.
        aligner.setMaxMemIntvOption(2000);
        aligner.setMaxSeedOccurencesOption(2000);
        // Other minor params to encourage more alignments to be found.
        aligner.setBandwidthOption(mTargetProbeLength);
        aligner.setSplitFactorOption(0.5f);
        aligner.setDropRatioOption(0.1f);
    }

    private static void logBwaMemParams(BwaMemAligner aligner)
    {
        LOGGER.debug("BWA-MEM options:");
        LOGGER.debug("  MinSeedLength: {}", aligner.getMinSeedLengthOption());
        LOGGER.debug("  SplitFactor: {}", aligner.getSplitFactorOption());
        LOGGER.debug("  SplitWidth: {}", aligner.getSplitWidthOption());
        LOGGER.debug("  MaxSeedOccurrences: {}", aligner.getMaxSeedOccurencesOption());
        LOGGER.debug("  MaxMemOccurrences: {}", aligner.getMaxMemIntvOption());
        LOGGER.debug("  DropRatio: {}", aligner.getDropRatioOption());
        LOGGER.debug("  Match: {}", aligner.getMatchScoreOption());
        LOGGER.debug("  Mismatch: {}", aligner.getMismatchPenaltyOption());
        LOGGER.debug("  IGapOpen: {}", aligner.getIGapOpenPenaltyOption());
        LOGGER.debug("  IGapExtend: {}", aligner.getIGapExtendPenaltyOption());
        LOGGER.debug("  DGapOpen: {}", aligner.getDGapOpenPenaltyOption());
        LOGGER.debug("  DGapExtend: {}", aligner.getDGapExtendPenaltyOption());
        LOGGER.debug("  Clip3: {}", aligner.getClip3PenaltyOption());
        LOGGER.debug("  Clip5: {}", aligner.getClip5PenaltyOption());
        LOGGER.debug("  Bandwidth: {}", aligner.getBandwidthOption());
        LOGGER.debug("  ZDrop: {}", aligner.getZDropOption());
        LOGGER.debug("  OutputScoreThreshold: {}", aligner.getOutputScoreThresholdOption());
        LOGGER.debug("  Flags: {}", aligner.getFlagOption());
        LOGGER.debug("  Threads: {}", aligner.getNThreadsOption());
    }
}
