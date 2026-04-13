package com.hartwig.hmftools.common.bwa;

public record BwaMemAlignParams(
        // Rewards and penalties are specified as positive numbers.
        int matchReward,    // -A
        int mismatchPenalty,    // -B
        int gapOpenPenalty,     // -O
        int gapExtendPenalty,   // -E
        int clipPenalty,    // -L
        int seedLengthMin,  // -k
        int seed3MaxOccurrence,     // -y and also known as max_mem_intv
        int memMaxOccurrence,   // -c
        float memReseedFactor,      // -r
        float chainOverlapFactor,   // -D
        int bandWidth,  // -w
        int zDropoff    // -d
)
{
    public BwaMemAlignParams
    {
        if(matchReward < 1 || mismatchPenalty < 0 || gapOpenPenalty < 0 || gapExtendPenalty < 0 || clipPenalty < 0)
        {
            throw new IllegalArgumentException("Invalid rewards/penalties");
        }
        if(seedLengthMin < 5)
        {
            // Not sure what the actual limit is, but you probably shouldn't try to go lower than this.
            throw new IllegalArgumentException("Invalid seedLengthMin");
        }
        if(seed3MaxOccurrence <= 0 || memMaxOccurrence <= 0 || memReseedFactor <= 0 || chainOverlapFactor <= 0 || bandWidth <= 0
                || zDropoff <= 0)
        {
            throw new IllegalArgumentException("Invalid parameters");
        }
    }

    public static final BwaMemAlignParams DEFAULT = new BwaMemAlignParams(1, 4, 6, 1, 5, 19, 20, 500, 1.5f, 0.5f, 100, 100);
}
