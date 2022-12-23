package com.hartwig.hmftools.common.isofox;

import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.rna.ImmutableGeneExpression;
import com.hartwig.hmftools.common.rna.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.common.rna.ImmutableRnaFusion;
import com.hartwig.hmftools.common.rna.ImmutableRnaStatistics;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class IsofoxTestFactory {

    private IsofoxTestFactory() {
    }

    @NotNull
    public static IsofoxData createMinimalIsofoxTestData() {
        return ImmutableIsofoxData.builder().summary(IsofoxTestFactory.rnaStatisticsBuilder().build()).build();
    }
    
    @NotNull
    public static ImmutableRnaStatistics.Builder rnaStatisticsBuilder() {
        return ImmutableRnaStatistics.builder()
                .totalFragments(0)
                .duplicateFragments(0)
                .splicedFragmentPerc(0D)
                .unsplicedFragmentPerc(0D)
                .altFragmentPerc(0D)
                .chimericFragmentPerc(0D)
                .readLength(0)
                .fragmentLength5thPercent(0D)
                .fragmentLength50thPercent(0D)
                .fragmentLength95thPercent(0D)
                .enrichedGenePercent(0D)
                .medianGCRatio(0D);
    }

    @NotNull
    public static ImmutableGeneExpression.Builder geneExpressionBuilder() {
        return ImmutableGeneExpression.builder()
                .geneName(Strings.EMPTY)
                .tpm(0D)
                .splicedFragments(0)
                .unsplicedFragments(0)
                .medianTpmCancer(0D)
                .percentileCancer(0D)
                .medianTpmCohort(0D)
                .percentileCohort(0D);
    }

    @NotNull
    public static ImmutableRnaFusion.Builder rnaFusionBuilder() {
        return ImmutableRnaFusion.builder()
                .name(Strings.EMPTY)
                .chromosomeUp(Strings.EMPTY)
                .chromosomeDown(Strings.EMPTY)
                .positionUp(0)
                .positionDown(0)
                .orientationUp((byte) 1)
                .orientationDown((byte) 0)
                .junctionTypeUp(Strings.EMPTY)
                .junctionTypeDown(Strings.EMPTY)
                .svType(StructuralVariantType.BND)
                .splitFragments(0)
                .realignedFrags(0)
                .discordantFrags(0)
                .depthUp(0)
                .depthDown(0)
                .maxAnchorLengthUp(0)
                .maxAnchorLengthDown(0)
                .cohortFrequency(0);
    }

    @NotNull
    public static ImmutableNovelSpliceJunction.Builder novelSpliceJunctionBuilder() {
        return ImmutableNovelSpliceJunction.builder()
                .chromosome(Strings.EMPTY)
                .junctionStart(0)
                .junctionEnd(0)
                .type(AltSpliceJunctionType.UNKNOWN)
                .fragmentCount(0)
                .depthStart(0)
                .depthEnd(0)
                .regionStart(AltSpliceJunctionContext.UNKNOWN)
                .regionEnd(AltSpliceJunctionContext.UNKNOWN)
                .basesStart(Strings.EMPTY)
                .basesEnd(Strings.EMPTY)
                .cohortFrequency(0);
    }
}
