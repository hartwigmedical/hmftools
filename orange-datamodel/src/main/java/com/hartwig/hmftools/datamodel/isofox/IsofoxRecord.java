package com.hartwig.hmftools.datamodel.isofox;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface IsofoxRecord
{
    @NotNull
    RnaStatistics summary();

    @NotNull
    List<GeneExpression> highExpressionGenes();

    @NotNull
    List<GeneExpression> lowExpressionGenes();

    @NotNull
    List<RnaFusion> fusions();

    @NotNull
    List<NovelSpliceJunction> novelSpliceJunctions();
}
