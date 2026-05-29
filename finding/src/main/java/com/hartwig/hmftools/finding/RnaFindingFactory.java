package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.RNA_REQUIRED;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.RnaStatistics;
import com.hartwig.hmftools.finding.datamodel.NovelSpliceJunction;
import com.hartwig.hmftools.finding.datamodel.NovelSpliceJunctionBuilder;
import com.hartwig.hmftools.finding.datamodel.RnaFusion;
import com.hartwig.hmftools.finding.datamodel.RnaFusionBuilder;
import com.hartwig.hmftools.finding.datamodel.RnaGeneExpression;
import com.hartwig.hmftools.finding.datamodel.RnaGeneExpressionBuilder;
import com.hartwig.hmftools.finding.datamodel.RnaStatisticsBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.util.FindingUtil;

import org.jetbrains.annotations.Nullable;

final class RnaFindingFactory
{
    static FindingItem<com.hartwig.hmftools.finding.datamodel.RnaStatistics> createRnaStatistics(
            @Nullable IsofoxRecord isofox, FindingStatus findingStatus)
    {
        if(isofox == null)
        {
            return FindingUtil.notAvailableFindingItem(Set.of(RNA_REQUIRED));
        }

        return FindingItemBuilder.<com.hartwig.hmftools.finding.datamodel.RnaStatistics>builder()
                .status(FindingUtil.somaticStatus(findingStatus))
                .finding(convertRnaStatistics(isofox.summary()))
                .build();
    }

    static FindingList<RnaGeneExpression> createRnaGeneExpressionFindings(@Nullable IsofoxRecord isofox,
            boolean highExpression, FindingStatus findingStatus)
    {
        if(isofox == null)
        {
            return rnaNotAvailableFindingList();
        }

        List<GeneExpression> geneExpressions = highExpression ? isofox.highExpressionGenes() : isofox.lowExpressionGenes();
        String expressionType = highExpression ? "HIGH" : "LOW";
        return FindingListBuilder.<RnaGeneExpression>builder()
                .status(FindingUtil.somaticStatus(findingStatus))
                .findings(geneExpressions.stream()
                        .map(geneExpression -> convertRnaGeneExpression(expressionType, geneExpression))
                        .sorted(RnaGeneExpression.COMPARATOR)
                        .toList())
                .build();
    }

    static FindingList<RnaFusion> createRnaFusionFindings(@Nullable IsofoxRecord isofox, FindingStatus findingStatus)
    {
        if(isofox == null)
        {
            return rnaNotAvailableFindingList();
        }

        return FindingListBuilder.<RnaFusion>builder()
                .status(FindingUtil.somaticStatus(findingStatus))
                .findings(isofox.fusions().stream()
                        .map(RnaFindingFactory::convertRnaFusion)
                        .sorted(RnaFusion.COMPARATOR)
                        .toList())
                .build();
    }

    static FindingList<NovelSpliceJunction> createNovelSpliceJunctionFindings(@Nullable IsofoxRecord isofox,
            FindingStatus findingStatus)
    {
        if(isofox == null)
        {
            return rnaNotAvailableFindingList();
        }

        return FindingListBuilder.<NovelSpliceJunction>builder()
                .status(FindingUtil.somaticStatus(findingStatus))
                .findings(isofox.novelSpliceJunctions().stream()
                        .map(RnaFindingFactory::convertNovelSpliceJunction)
                        .sorted(NovelSpliceJunction.COMPARATOR)
                        .toList())
                .build();
    }

    private static com.hartwig.hmftools.finding.datamodel.RnaStatistics convertRnaStatistics(RnaStatistics statistics)
    {
        return RnaStatisticsBuilder.builder()
                .findingKey(FindingKeys.rnaStatistics())
                .qcStatus(statistics.qcStatus().stream()
                        .map(status -> com.hartwig.hmftools.finding.datamodel.RnaStatistics.QcStatus.valueOf(status.name()))
                        .collect(Collectors.toSet()))
                .totalFragments(statistics.totalFragments())
                .duplicateFragments(statistics.duplicateFragments())
                .splicedFragmentPercent(statistics.splicedFragmentPerc())
                .unsplicedFragmentPercent(statistics.unsplicedFragmentPerc())
                .altFragmentPercent(statistics.altFragmentPerc())
                .chimericFragmentPercent(statistics.chimericFragmentPerc())
                .build();
    }

    private static RnaGeneExpression convertRnaGeneExpression(String expressionType, GeneExpression geneExpression)
    {
        return RnaGeneExpressionBuilder.builder()
                .findingKey(FindingKeys.rnaGeneExpression(expressionType, geneExpression))
                .gene(geneExpression.gene())
                .tpm(geneExpression.tpm())
                .medianTpmCohort(geneExpression.medianTpmCohort())
                .percentileCohort(geneExpression.percentileCohort())
                .medianTpmCancer(geneExpression.medianTpmCancer())
                .percentileCancer(geneExpression.percentileCancer())
                .build();
    }

    private static RnaFusion convertRnaFusion(com.hartwig.hmftools.datamodel.isofox.RnaFusion fusion)
    {
        return RnaFusionBuilder.builder()
                .findingKey(FindingKeys.rnaFusion(fusion))
                .geneStart(fusion.geneStart())
                .geneEnd(fusion.geneEnd())
                .chromosomeStart(fusion.chromosomeStart())
                .chromosomeEnd(fusion.chromosomeEnd())
                .positionStart(fusion.positionStart())
                .positionEnd(fusion.positionEnd())
                .junctionTypeStart(fusion.junctionTypeStart())
                .junctionTypeEnd(fusion.junctionTypeEnd())
                .knownType(RnaFusion.KnownType.valueOf(fusion.knownType().name()))
                .structuralVariantType(RnaFusion.StructuralVariantType.valueOf(fusion.svType().name()))
                .splitFragments(fusion.splitFragments())
                .realignedFragments(fusion.realignedFrags())
                .discordantFragments(fusion.discordantFrags())
                .depthStart(fusion.depthStart())
                .depthEnd(fusion.depthEnd())
                .cohortFrequency(fusion.cohortFrequency())
                .build();
    }

    private static NovelSpliceJunction convertNovelSpliceJunction(
            com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction spliceJunction)
    {
        return NovelSpliceJunctionBuilder.builder()
                .findingKey(FindingKeys.novelSpliceJunction(spliceJunction))
                .gene(spliceJunction.gene())
                .chromosome(spliceJunction.chromosome())
                .junctionStart(spliceJunction.junctionStart())
                .junctionEnd(spliceJunction.junctionEnd())
                .type(NovelSpliceJunction.Type.valueOf(spliceJunction.type().name()))
                .exonStart(spliceJunction.exonStart())
                .exonEnd(spliceJunction.exonEnd())
                .fragmentCount(spliceJunction.fragmentCount())
                .depthStart(spliceJunction.depthStart())
                .depthEnd(spliceJunction.depthEnd())
                .regionStart(NovelSpliceJunction.Context.valueOf(spliceJunction.regionStart().name()))
                .regionEnd(NovelSpliceJunction.Context.valueOf(spliceJunction.regionEnd().name()))
                .cohortFrequency(spliceJunction.cohortFrequency())
                .build();
    }

    private static <T> FindingList<T> rnaNotAvailableFindingList()
    {
        return FindingListBuilder.<T>builder()
                .status(FindingUtil.notAvailableStatus(Set.of(RNA_REQUIRED)))
                .findings(List.of())
                .build();
    }
}
