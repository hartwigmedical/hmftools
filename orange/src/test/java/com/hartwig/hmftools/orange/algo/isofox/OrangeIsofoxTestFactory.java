package com.hartwig.hmftools.orange.algo.isofox;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableGeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.ImmutableRnaFusion;
import com.hartwig.hmftools.datamodel.isofox.ImmutableRnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.RnaStatistics;
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
    public static ImmutableRnaStatistics.Builder rnaStatisticsBuilder()
    {
        com.hartwig.hmftools.common.rna.RnaStatistics common = IsofoxTestFactory.rnaStatisticsBuilder().build();
        RnaStatistics converted = IsofoxConversion.convert(common);
        return ImmutableRnaStatistics.builder().from(converted);
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
