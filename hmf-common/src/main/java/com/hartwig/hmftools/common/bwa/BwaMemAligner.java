package com.hartwig.hmftools.common.bwa;

import static org.broadinstitute.hellbender.utils.bwa.BwaMemAligner.MEM_F_ALL;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

public class BwaMemAligner implements IBwaMemAligner
{
    private final org.broadinstitute.hellbender.utils.bwa.BwaMemAligner mAligner;

    private static final Logger LOGGER = LogManager.getLogger(BwaMemAligner.class);

    public BwaMemAligner(final BwaMemIndex index, final BwaMemAlignParams alignParams, boolean allAlignments, int minAlignScore,
            int threads)
    {
        if(threads < 1)
        {
            throw new IllegalArgumentException("Invalid threads: " + threads);
        }
        if(minAlignScore < 0)
        {
            throw new IllegalArgumentException("Invalid minAlignScore: " + minAlignScore);
        }

        mAligner = new org.broadinstitute.hellbender.utils.bwa.BwaMemAligner(index);
        applyOptions(mAligner, alignParams, allAlignments, minAlignScore, threads);
        logOptions(mAligner);
    }

    public List<List<BwaMemAlignment>> alignSequences(List<byte[]> sequences)
    {
        return mAligner.alignSeqs(sequences);
    }

    private static void applyOptions(org.broadinstitute.hellbender.utils.bwa.BwaMemAligner aligner, BwaMemAlignParams alignParams,
            boolean allAlignments, int minAlignScore, int threads)
    {
        aligner.setMatchScoreOption(alignParams.matchReward());
        aligner.setMismatchPenaltyOption(alignParams.mismatchPenalty());
        aligner.setIGapOpenPenaltyOption(alignParams.gapOpenPenalty());
        aligner.setDGapOpenPenaltyOption(alignParams.gapOpenPenalty());
        aligner.setIGapExtendPenaltyOption(alignParams.gapExtendPenalty());
        aligner.setDGapExtendPenaltyOption(alignParams.gapExtendPenalty());
        aligner.setClip3PenaltyOption(alignParams.clipPenalty());
        aligner.setClip5PenaltyOption(alignParams.clipPenalty());
        // I'm not sure if this is strictly necessary, since the other params should set the scoring.
        // But it was done is ESVEE, so to ensure backwards compatibility, we will do it.
        updateScoringMatrix(aligner);

        aligner.setMinSeedLengthOption(alignParams.seedLengthMin());
        aligner.setMaxMemIntvOption(alignParams.seed3MaxOccurrence());
        aligner.setMaxSeedOccurencesOption(alignParams.memMaxOccurrence());
        aligner.setSplitFactorOption(alignParams.memReseedFactor());
        aligner.setDropRatioOption(alignParams.chainOverlapFactor());
        aligner.setBandwidthOption(alignParams.bandWidth());
        aligner.setZDropOption(alignParams.zDropoff());

        if(allAlignments)
        {
            aligner.setFlagOption(aligner.getFlagOption() | MEM_F_ALL);
        }
        else
        {
            aligner.setFlagOption(aligner.getFlagOption() & ~MEM_F_ALL);
        }
        aligner.setOutputScoreThresholdOption(minAlignScore);

        aligner.setNThreadsOption(threads);
    }

    private static void updateScoringMatrix(org.broadinstitute.hellbender.utils.bwa.BwaMemAligner aligner)
    {
        // TODO: what is this actually achieving?

        int matrixSize = 5;
        int matchScore = aligner.getMatchScoreOption();
        int mismatchPenalty = aligner.getMismatchPenaltyOption();
        byte[] scoringMatrix = new byte[matrixSize * matrixSize];
        int k = 0;

        for(int i = 0; i < matrixSize - 1; i++)
        {
            for(int j = 0; j < matrixSize - 1; j++)
            {
                scoringMatrix[k++] = (byte) (i == j ? matchScore : -mismatchPenalty);
            }
            scoringMatrix[k++] = -1;
        }

        for(int j = 0; j < matrixSize; j++)
        {
            scoringMatrix[k++] = -1;
        }

        aligner.setScoringMatrixOption(scoringMatrix);
    }

    private static void logOptions(org.broadinstitute.hellbender.utils.bwa.BwaMemAligner aligner)
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
