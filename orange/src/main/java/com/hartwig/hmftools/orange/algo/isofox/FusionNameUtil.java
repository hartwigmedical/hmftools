package com.hartwig.hmftools.orange.algo.isofox;

import com.hartwig.hmftools.common.rna.RnaFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FusionNameUtil
{
    private static final String RNA_FUSION_NAME_DELIMITER = "_";

    @Nullable
    public static String geneUp(@NotNull RnaFusion rnaFusion)
    {
        int split = rnaFusion.name().indexOf(RNA_FUSION_NAME_DELIMITER);
        return split > 0 ? rnaFusion.name().substring(0, split) : null;
    }

    @Nullable
    public static String geneDown(@NotNull RnaFusion rnaFusion)
    {
        int split = rnaFusion.name().indexOf(RNA_FUSION_NAME_DELIMITER);
        return split >= 0 && split < rnaFusion.name().length() - 1 ? rnaFusion.name().substring(split + 1) : null;
    }
}
