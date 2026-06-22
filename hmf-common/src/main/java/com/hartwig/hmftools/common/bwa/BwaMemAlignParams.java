package com.hartwig.hmftools.common.bwa;

public record BwaMemAlignParams(
        // Rewards and penalties are specified as positive numbers.

        // Reward for matching bases, except for N. CLI: -A
        int matchReward,
        // Penalty for mismatched bases, except if one is N. CLI: -B
        int mismatchPenalty,
        // Penalty for when either base is N. No CLI equivalent.
        int uncertainBasePenalty,
        // Penalty for opening an insert or delete. CLI: -O
        int gapOpenPenalty,
        // Penalty for adding a base to an insert or delete. CLI: -E
        int gapExtendPenalty,
        // Penalty for terminating the alignment early. CLI: -L
        int clipPenalty,
        // Length of exact match seed. CLI: -k
        int seedLengthMin,
        // In the 3rd round of seeding, seeds with more than this many occurrences are dropped. Also known as max_mem_intv. CLI: -y
        int seed3MaxOccurrence,
        // CLI: -c
        int memMaxOccurrence,
        // CLI: -r
        float memReseedFactor,
        // CLI: -D
        float chainOverlapFactor,
        // CLI: -w
        int bandWidth,
        // CLI: -d
        int zDropoff,
        // Minimum alignment score. CLI: -T
        int minAlignScore
)
{
    public BwaMemAlignParams
    {
        if(matchReward < 1 || mismatchPenalty < 0 || uncertainBasePenalty < 0 || gapOpenPenalty < 0 || gapExtendPenalty < 0
                || clipPenalty < 0)
        {
            throw new IllegalArgumentException("Invalid rewards/penalties");
        }

        if(seedLengthMin < 5)
        {
            // Not sure what the actual limit is, but you probably shouldn't try to go lower than this.
            throw new IllegalArgumentException("Invalid seedLengthMin");
        }

        if(seed3MaxOccurrence <= 0 || memMaxOccurrence <= 0 || memReseedFactor <= 0 || chainOverlapFactor <= 0 || bandWidth <= 0
                || zDropoff <= 0 || minAlignScore < 0)
        {
            throw new IllegalArgumentException("Invalid parameters");
        }
    }

    public BwaMemAlignParams withMismatchPenalty(int value)
    {
        return new BwaMemAlignParams(matchReward, value, uncertainBasePenalty, gapOpenPenalty, gapExtendPenalty, clipPenalty,
                seedLengthMin, seed3MaxOccurrence, memMaxOccurrence, memReseedFactor, chainOverlapFactor, bandWidth, zDropoff, minAlignScore);
    }

    public BwaMemAlignParams withGapOpenPenalty(int value)
    {
        return new BwaMemAlignParams(matchReward, mismatchPenalty, uncertainBasePenalty, value, gapExtendPenalty, clipPenalty,
                seedLengthMin, seed3MaxOccurrence, memMaxOccurrence, memReseedFactor, chainOverlapFactor, bandWidth, zDropoff, minAlignScore);
    }

    public BwaMemAlignParams withMinAlignScore(int value)
    {
        return new BwaMemAlignParams(matchReward, mismatchPenalty, uncertainBasePenalty, gapOpenPenalty, gapExtendPenalty, clipPenalty,
                seedLengthMin, seed3MaxOccurrence, memMaxOccurrence, memReseedFactor, chainOverlapFactor, bandWidth, zDropoff, value);
    }

    public static final BwaMemAlignParams DEFAULT = new BwaMemAlignParams(
            1, 4, 1, 6, 1, 5, 19,
            20, 500, 1.5f, 0.5f, 100, 100, 30);
}
