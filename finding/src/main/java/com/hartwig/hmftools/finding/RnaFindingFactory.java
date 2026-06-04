package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.RNA_REQUIRED;
import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.RNA_SAMPLE_QUALITY_CONTROL;

import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.RnaQCStatus;
import com.hartwig.hmftools.finding.datamodel.NovelSpliceJunction;
import com.hartwig.hmftools.finding.datamodel.NovelSpliceJunctionBuilder;
import com.hartwig.hmftools.finding.datamodel.RnaFusion;
import com.hartwig.hmftools.finding.datamodel.RnaFusionBuilder;
import com.hartwig.hmftools.finding.datamodel.RnaGeneExpression;
import com.hartwig.hmftools.finding.datamodel.RnaGeneExpressionBuilder;
import com.hartwig.hmftools.finding.datamodel.RnaStatistics;
import com.hartwig.hmftools.finding.datamodel.RnaStatisticsBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;
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
                .status(rnaStatus(isofox, findingStatus))
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
                .status(rnaStatus(isofox, findingStatus))
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
                .status(rnaStatus(isofox, findingStatus))
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
                .status(rnaStatus(isofox, findingStatus))
                .findings(isofox.novelSpliceJunctions().stream()
                        .map(RnaFindingFactory::convertNovelSpliceJunction)
                        .sorted(NovelSpliceJunction.COMPARATOR)
                        .toList())
                .build();
    }

    private static RnaStatistics convertRnaStatistics(com.hartwig.hmftools.datamodel.isofox.RnaStatistics statistics)
    {
        SortedSet<RnaStatistics.QcStatus> qcErrors = new TreeSet<>();
        SortedSet<RnaStatistics.QcStatus> qcWarnings = new TreeSet<>();

        statistics.qcStatus().forEach(qcStatus -> {
            switch(qcStatus)
            {
                case PASS:
                    break;
                case FAIL_LOW_COVERAGE:
                    qcErrors.add(RnaStatistics.QcStatus.FAIL_LOW_COVERAGE);
                    break;
                case WARN_DUPLICATE_RATE:
                    qcWarnings.add(RnaStatistics.QcStatus.WARN_DUPLICATE_RATE);
                    break;
                case WARN_SPLICED_GENE_COVERAGE:
                    qcWarnings.add(RnaStatistics.QcStatus.WARN_SPLICED_GENE_COVERAGE);
                    break;
            }
        });

        return RnaStatisticsBuilder.builder()
                .findingKey(FindingKeys.rnaStatistics())
                .errors(qcErrors)
                .warnings(qcWarnings)
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
                .geneUp(fusion.geneStart())
                .geneDown(fusion.geneEnd())
                .chromosomeUp(fusion.chromosomeStart())
                .chromosomeDown(fusion.chromosomeEnd())
                .positionUp(fusion.positionStart())
                .positionDown(fusion.positionEnd())
                .junctionTypeUp(fusion.junctionTypeStart())
                .junctionTypeDown(fusion.junctionTypeEnd())
                .knownType(RnaFusion.KnownType.valueOf(fusion.knownType().name()))
                .structuralVariantType(RnaFusion.StructuralVariantType.valueOf(fusion.svType().name()))
                .splitFragments(fusion.splitFragments())
                .realignedFragments(fusion.realignedFrags())
                .discordantFragments(fusion.discordantFrags())
                .depthUp(fusion.depthStart())
                .depthDown(fusion.depthEnd())
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

    private static FindingStatus rnaStatus(IsofoxRecord isofox, FindingStatus findingStatus)
    {
        FindingStatus somaticStatus = FindingUtil.somaticStatus(findingStatus);
        if(isofox.summary().qcStatus().stream().anyMatch(RnaFindingFactory::isFailingQcStatus))
        {
            return FindingStatusBuilder.builder(somaticStatus)
                    .status(FindingStatus.Status.NOT_RELIABLE)
                    .errors(FindingUtil.addIssues(somaticStatus.errors(), Set.of(RNA_SAMPLE_QUALITY_CONTROL)))
                    .build();
        }
        return somaticStatus;
    }

    private static boolean isFailingQcStatus(RnaQCStatus qcStatus)
    {
        return qcStatus.name().startsWith("FAIL");
    }
}
