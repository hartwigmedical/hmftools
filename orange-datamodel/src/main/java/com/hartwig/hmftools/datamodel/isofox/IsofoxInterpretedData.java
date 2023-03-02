package com.hartwig.hmftools.datamodel.isofox;

import com.hartwig.hmftools.datamodel.rna.GeneExpression;
import com.hartwig.hmftools.datamodel.rna.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.rna.RnaFusion;
import com.hartwig.hmftools.datamodel.rna.RnaStatistics;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class IsofoxInterpretedData {

    @NotNull
    public abstract RnaStatistics summary();

    @NotNull
    public abstract List<GeneExpression> allGeneExpressions();

    @NotNull
    public abstract List<GeneExpression> reportableHighExpression();

    @NotNull
    public abstract List<GeneExpression> reportableLowExpression();

    @NotNull
    public abstract List<RnaFusion> allFusions();

    @NotNull
    public abstract List<RnaFusion> reportableNovelKnownFusions();

    @NotNull
    public abstract List<RnaFusion> reportableNovelPromiscuousFusions();

    @NotNull
    public abstract List<NovelSpliceJunction> allNovelSpliceJunctions();

    @NotNull
    public abstract List<NovelSpliceJunction> reportableSkippedExons();

    @NotNull
    public abstract List<NovelSpliceJunction> reportableNovelExonsIntrons();

}
