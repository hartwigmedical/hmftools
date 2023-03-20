package com.hartwig.hmftools.datamodel.isofox;

import com.hartwig.hmftools.datamodel.rna.GeneExpression;
import com.hartwig.hmftools.datamodel.rna.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.rna.RnaFusion;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Value.Immutable
@Value.Style(passAnnotations = {NotNull.class, Nullable.class})
public interface IsofoxInterpretedData {

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
