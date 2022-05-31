package com.hartwig.hmftools.orange.isofox;

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
public interface IsofoxInterpretedData {

    @NotNull
    RnaStatistics summary();

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
