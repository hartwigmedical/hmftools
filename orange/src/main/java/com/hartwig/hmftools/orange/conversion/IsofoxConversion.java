package com.hartwig.hmftools.orange.conversion;

import static com.hartwig.hmftools.orange.algo.isofox.FusionNameUtil.geneDown;
import static com.hartwig.hmftools.orange.algo.isofox.FusionNameUtil.geneUp;

import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaQcFilter;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionContext;
import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionType;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableGeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.ImmutableRnaFusion;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.RnaQCStatus;
import com.hartwig.hmftools.datamodel.isofox.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

public final class IsofoxConversion
{
    @NotNull
    public static IsofoxRnaStatistics convert(@NotNull RnaStatistics rnaStatistics)
    {
        return ImmutableIsofoxRnaStatistics.builder()
                .totalFragments(rnaStatistics.totalFragments())
                .duplicateFragments(rnaStatistics.duplicateFragments())
                .qcStatus(ConversionUtil.mapToIterable(rnaStatistics.qcStatus(), IsofoxConversion::convert))
                .build();
    }

    @NotNull
    public static RnaQCStatus convert(@NotNull RnaQcFilter filter)
    {
        return RnaQCStatus.valueOf(filter.name());
    }

    @NotNull
    public static GeneExpression convert(@NotNull com.hartwig.hmftools.common.rna.GeneExpression geneExpression)
    {
        return ImmutableGeneExpression.builder()
                .gene(geneExpression.geneName())
                .tpm(geneExpression.tpm())
                .medianTpmCancer(geneExpression.medianTpmCancer() >= 0 ? geneExpression.medianTpmCancer() : null)
                .percentileCancer(geneExpression.percentileCancer() >= 0 ? geneExpression.percentileCancer() : null)
                .medianTpmCohort(geneExpression.medianTpmCohort())
                .percentileCohort(geneExpression.percentileCohort())
                .build();
    }

    @NotNull
    public static com.hartwig.hmftools.datamodel.isofox.RnaFusion convert(@NotNull RnaFusion rnaFusion)
    {
        return ImmutableRnaFusion.builder()
                .geneStart(geneUp(rnaFusion))
                .geneEnd(geneDown(rnaFusion))
                .chromosomeStart(rnaFusion.chromosomeUp())
                .chromosomeEnd(rnaFusion.chromosomeDown())
                .positionStart(rnaFusion.positionUp())
                .positionEnd(rnaFusion.positionDown())
                .junctionTypeStart(rnaFusion.junctionTypeUp())
                .junctionTypeEnd(rnaFusion.junctionTypeDown())
                .svType(StructuralVariantType.valueOf(rnaFusion.svType().name()))
                .splitFragments(rnaFusion.splitFragments())
                .realignedFrags(rnaFusion.realignedFrags())
                .discordantFrags(rnaFusion.discordantFrags())
                .depthStart(rnaFusion.depthUp())
                .depthEnd(rnaFusion.depthDown())
                .cohortFrequency(rnaFusion.cohortFrequency())
                .build();
    }

    @NotNull
    public static com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction convert(@NotNull NovelSpliceJunction novelSpliceJunction)
    {
        return ImmutableNovelSpliceJunction.builder()
                .gene(novelSpliceJunction.geneName())
                .chromosome(novelSpliceJunction.chromosome())
                .junctionStart(novelSpliceJunction.junctionStart())
                .junctionEnd(novelSpliceJunction.junctionEnd())
                .type(AltSpliceJunctionType.valueOf(novelSpliceJunction.type().name()))
                .fragmentCount(novelSpliceJunction.fragmentCount())
                .depthStart(novelSpliceJunction.depthStart())
                .depthEnd(novelSpliceJunction.depthEnd())
                .regionStart(AltSpliceJunctionContext.valueOf(novelSpliceJunction.regionStart().name()))
                .regionEnd(AltSpliceJunctionContext.valueOf(novelSpliceJunction.regionEnd().name()))
                .cohortFrequency(novelSpliceJunction.cohortFrequency())
                .build();
    }
}
