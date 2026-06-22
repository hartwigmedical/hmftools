package com.hartwig.hmftools.common.bwa;

import org.jetbrains.annotations.Nullable;

public record BwaMemAlignerConfig(
        String indexPath,
        BwaMemAlignParams alignParams,
        boolean allAlignments,
        int threads,
        @Nullable Integer batchSize
)
{
    public BwaMemAlignerConfig
    {
        if(threads < 1)
        {
            throw new IllegalArgumentException("Invalid threads: " + threads);
        }
        if(batchSize != null && batchSize < 1)
        {
            throw new IllegalArgumentException("Invalid batchSize: " + batchSize);
        }
    }

    public BwaMemAlignerConfig(final String indexPath, BwaMemAlignParams align)
    {
        this(indexPath, align, false, 1, null);
    }

    public BwaMemAlignerConfig(final String indexPath)
    {
        this(indexPath, BwaMemAlignParams.DEFAULT);
    }
}
