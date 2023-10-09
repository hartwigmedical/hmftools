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
    IsofoxRnaStatistics summary();

    @NotNull
    List<GeneExpression> allGeneExpressions();

    @NotNull
    List<GeneExpression> reportableHighExpression();

    @NotNull
    List<GeneExpression> reportableLowExpression();

    @NotNull
    List<RnaFusion> allFusions();

    @NotNull
    List<RnaFusion> reportableNovelKnownFusions();

    @NotNull
    List<RnaFusion> reportableNovelPromiscuousFusions();

    @NotNull
    List<NovelSpliceJunction> allNovelSpliceJunctions();

    @NotNull
    List<NovelSpliceJunction> reportableSkippedExons();

    @NotNull
    List<NovelSpliceJunction> reportableNovelExonsIntrons();

}
