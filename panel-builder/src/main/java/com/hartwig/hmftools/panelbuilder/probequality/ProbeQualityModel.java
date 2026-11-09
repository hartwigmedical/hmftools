package com.hartwig.hmftools.panelbuilder.probequality;

import static java.lang.Math.min;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.bwa.BwaMemAlignParams;
import com.hartwig.hmftools.common.bwa.BwaMemAlignerConfig;
import com.hartwig.hmftools.common.bwa.IBwaMemAligner;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
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
    // BWA-MEM per-base match reward; a full-length perfect alignment scores mTargetProbeLength * mMatchScore.
    private final int mMatchScore;
    // Maps a BWA alignment's reference contig index (getRefId) to its chromosome name, so an alignment can be tested against a probe's source
    // regions. Built from the reference sequence dictionary, whose order matches the BWA index.
    private final Map<Integer, String> mRefIdToChromosome;

    // Use in production; the constructor is for injecting a test aligner.
    public static ProbeQualityModel create(
            final String bwaIndexImage, final int threads, final int targetProbeLength, final int matchScoreThreshold,
            final int matchScoreOffset, final Map<Integer, String> refIdToChromosome)
    {
        BwaMemAlignerConfig config = new BwaMemAlignerConfig(
                bwaIndexImage, alignParams(targetProbeLength, matchScoreThreshold), true, threads, null);
        return new ProbeQualityModel(
                new BwaMemAligner(config), targetProbeLength, matchScoreThreshold, matchScoreOffset, refIdToChromosome);
    }

    public ProbeQualityModel(final IBwaMemAligner aligner, final int targetProbeLength, final int matchScoreThreshold,
            final int matchScoreOffset, final Map<Integer, String> refIdToChromosome)
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
        mMatchScore = BwaMemAlignParams.DEFAULT.matchReward();
        mRefIdToChromosome = refIdToChromosome;
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

    // Builds the refId -> chromosome map from a reference genome sequence dictionary, whose order matches the BWA index.
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

    // Compute probe qualities. probeSourceRegions[i] is the reference region(s) probe i is built from (its intended on-target capture loci):
    // one region for a normal probe, several for a constructed probe (variant/SV, RNA spliced). Alignments landing on those are excluded from
    // the off-target risk; the rest count as off-target, normalised against a theoretical full-length match.
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

        List<List<BwaMemAlignment>> alignments = mAligner.alignSequences(probes);
        return IntStream.range(0, alignments.size())
                .mapToObj(i -> computeFromAlignments(alignments.get(i), probeSourceRegions.get(i)))
                .toList();
    }

    private Result computeFromAlignments(final List<BwaMemAlignment> alignments, final List<ChrBaseRegion> sourceRegions)
    {
        // Normalise against a theoretical perfect full-length match. Alignments landing on the probe's source region(s) are the intended
        // on-target captures and are excluded; every other alignment above the threshold is off-target.
        int targetScore = mTargetProbeLength * mMatchScore;
        List<Integer> offTarget = alignments.stream()
                .filter(alignment -> !isOnTarget(alignment, sourceRegions))
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

    // Whether an alignment lands on one of the probe's source regions (an intended on-target capture), resolving the alignment's reference
    // contig index to a chromosome name via the sequence dictionary.
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
}
