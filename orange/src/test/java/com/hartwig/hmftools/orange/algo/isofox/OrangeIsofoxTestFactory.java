package com.hartwig.hmftools.orange.algo.isofox;

import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableGeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.ImmutableRnaFusion;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.orange.conversion.IsofoxConversion;

import org.jetbrains.annotations.NotNull;

public class OrangeIsofoxTestFactory
{

    private OrangeIsofoxTestFactory()
    {
    }

    @NotNull
    public static ImmutableIsofoxRnaStatistics.Builder rnaStatisticsBuilder()
    {
        RnaStatistics common = IsofoxTestFactory.rnaStatisticsBuilder().build();
        IsofoxRnaStatistics converted = IsofoxConversion.convert(common);
        return ImmutableIsofoxRnaStatistics.builder().from(converted);
    }

    @NotNull
    public static ImmutableGeneExpression.Builder geneExpressionBuilder()
    {
        com.hartwig.hmftools.common.rna.GeneExpression common = IsofoxTestFactory.geneExpressionBuilder().build();
        GeneExpression converted = IsofoxConversion.convert(common);
        return ImmutableGeneExpression.builder().from(converted);
    }

    @NotNull
    public static ImmutableRnaFusion.Builder rnaFusionBuilder()
    {
        com.hartwig.hmftools.common.rna.RnaFusion common = IsofoxTestFactory.rnaFusionBuilder().name("START_END").build();
        RnaFusion converted = IsofoxConversion.convert(common);
        return ImmutableRnaFusion.builder().from(converted);
    }

    @NotNull
    public static ImmutableNovelSpliceJunction.Builder novelSpliceJunctionBuilder()
    {
        com.hartwig.hmftools.common.rna.NovelSpliceJunction common = IsofoxTestFactory.novelSpliceJunctionBuilder().build();
        NovelSpliceJunction converted = IsofoxConversion.convert(common);
        return ImmutableNovelSpliceJunction.builder().from(converted);
    }
}
