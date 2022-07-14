package com.hartwig.hmftools.common.isofox;

import java.util.List;

import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaStatistics;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface IsofoxData
{
    @NotNull
    RnaStatistics summary();

    @NotNull
    List<GeneExpression> geneExpressions();

    @NotNull
    List<RnaFusion> fusions();

    @NotNull
    List<NovelSpliceJunction> novelSpliceJunctions();

}
