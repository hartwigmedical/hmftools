package com.hartwig.hmftools.common.bwa;

import static com.hartwig.hmftools.common.bwa.BwaUtils.LIBBWA_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.common.utils.Streams.partitionStream;

import static org.broadinstitute.hellbender.utils.bwa.BwaMemAligner.MEM_F_ALL;

import java.util.List;
import java.util.stream.Stream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.jetbrains.annotations.Nullable;

public class BwaMemAligner implements IBwaMemAligner
{
    private final org.broadinstitute.hellbender.utils.bwa.BwaMemAligner mAligner;
    @Nullable
    private final Integer mBatchSize;

    private static final Logger LOGGER = LogManager.getLogger(BwaMemAligner.class);

    // Must be called once to load the BWA-MEM DLL before using any functionality.
    public static void initLibrary(@Nullable final String path)
    {
        // Probably a mistake to try to initialise it multiple times, even though technically it may work.
        if(System.getProperties().containsKey(LIBBWA_PATH))
        {
            throw new RuntimeException("BWA-MEM library already loaded");
        }
        loadAlignerLibrary(path);
    }

    public BwaMemAligner(final BwaMemIndex index, final BwaMemAlignParams alignParams, boolean allAlignments, int minAlignScore,
            int threads, @Nullable Integer batchSize)
    {
        if(minAlignScore < 0)
        {
            throw new IllegalArgumentException("Invalid minAlignScore: " + minAlignScore);
        }
        if(threads < 1)
        {
            throw new IllegalArgumentException("Invalid threads: " + threads);
        }
        if(batchSize != null && batchSize < 1)
        {
            throw new IllegalArgumentException("Invalid batchSize: " + batchSize);
        }

        mBatchSize = batchSize;

        LOGGER.debug("Creating BWA-MEM aligner");
        mAligner = new org.broadinstitute.hellbender.utils.bwa.BwaMemAligner(index);
        applyOptions(mAligner, alignParams, allAlignments, minAlignScore, threads);
        logOptions(mAligner);
    }

    public List<List<BwaMemAlignment>> alignSequences(List<byte[]> sequences)
    {
        return runBatchedAlignment(sequences);
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
        // Setting the match and mismatch options doesn't do anything because BWA uses an internal matrix of base-to-base scores.
        // Usually when called from the CLI, BWA will update this matrix accordingly. But via the JNI wrapper, that code won't trigger, so
        // we have to do it ourselves.
        setScoringMatrix(aligner, alignParams);

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

    private static void setScoringMatrix(org.broadinstitute.hellbender.utils.bwa.BwaMemAligner aligner, BwaMemAlignParams alignParams)
    {
        // Set the base-to-base scoring.
        // If the bases are equal, use the match reward.
        // If the bases are unequal, use the mismatch penalty.
        // If either base is uncertain (N), use the uncertain base penalty.

        byte matchScore = (byte) alignParams.matchReward();
        byte mismatchScore = (byte) -alignParams.mismatchPenalty();
        byte uncertainScore = (byte) -alignParams.uncertainBasePenalty();

        int matrixSize = 5; // 4 nucleotides + uncertain
        byte[] scoringMatrix = new byte[matrixSize * matrixSize];
        int k = 0;

        for(int i = 0; i < matrixSize - 1; i++)
        {
            for(int j = 0; j < matrixSize - 1; j++)
            {
                scoringMatrix[k++] = i == j ? matchScore : mismatchScore;
            }
            scoringMatrix[k++] = uncertainScore;
        }

        for(int j = 0; j < matrixSize; j++)
        {
            scoringMatrix[k++] = uncertainScore;
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

    private List<List<BwaMemAlignment>> runBatchedAlignment(List<byte[]> queries)
    {
        int batchSize = mBatchSize == null ? queries.size() : mBatchSize;
        while(true)
        {
            try
            {
                Stream<List<byte[]>> batches =
                        batchSize < queries.size() ? partitionStream(queries.stream(), batchSize) : Stream.of(queries);
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
                    LOGGER.warn("Aligning sequences with batch size {} failed. Trying again with batch size {}",
                            batchSize, newBatchSize);
                    batchSize = newBatchSize;
                    // Loop and try again with the new batch size.
                }
                else
                {
                    throw e;
                }
            }
        }
    }
}
