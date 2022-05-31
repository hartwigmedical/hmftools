package com.hartwig.hmftools.orange.algo.isofox;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.isofox.IsofoxData;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class IsofoxInterpreter {

    private static final String RNA_FUSION_NAME_DELIMITER = "_";

    private IsofoxInterpreter() {
    }

    @NotNull
    public static IsofoxInterpretedData interpret(@NotNull IsofoxData isofox, @NotNull LinxData linx, @NotNull List<DriverGene> driverGenes,
            @NotNull KnownFusionCache knownFusionCache) {
        return ImmutableIsofoxInterpretedData.builder()
                .summary(isofox.summary())
                .allGeneExpressions(isofox.geneExpressions())
                .reportableHighExpression(selectHighExpressionGenes(isofox.geneExpressions(), driverGenes))
                .reportableLowExpression(selectLowExpressionGenes(isofox.geneExpressions(), driverGenes))
                .allFusions(isofox.fusions())
                .reportableNovelKnownFusions(selectNovelKnownFusions(isofox.fusions(), linx.reportableFusions(), knownFusionCache))
                .reportableNovelPromiscuousFusions(selectNovelPromiscuousFusions(isofox.fusions(),
                        linx.reportableFusions(),
                        knownFusionCache))
                .allNovelSpliceJunctions(isofox.novelSpliceJunctions())
                .reportableSkippedExons(selectSkippedExons(isofox.novelSpliceJunctions(), linx.reportableFusions(), knownFusionCache))
                .reportableNovelExonsIntrons(selectNovelExonsIntrons(isofox.novelSpliceJunctions(), driverGenes))
                .build();
    }

    @NotNull
    private static List<GeneExpression> selectHighExpressionGenes(@NotNull List<GeneExpression> expressions,
            @NotNull List<DriverGene> driverGenes) {
        Set<String> oncogenes = extractGenesOfType(driverGenes, DriverCategory.ONCO);

        List<GeneExpression> result = Lists.newArrayList();
        for (GeneExpression expression : expressions) {
            if (oncogenes.contains(expression.geneName()) && expression.percentileCohort() > 0.9 && expression.percentileCancer() > 0.9) {
                result.add(expression);
            }
        }

        return result;
    }

    @NotNull
    private static List<GeneExpression> selectLowExpressionGenes(@NotNull List<GeneExpression> expressions,
            @NotNull List<DriverGene> driverGenes) {
        Set<String> suppressors = extractGenesOfType(driverGenes, DriverCategory.TSG);

        List<GeneExpression> result = Lists.newArrayList();
        for (GeneExpression expression : expressions) {
            if (suppressors.contains(expression.geneName()) && expression.percentileCohort() < 0.05
                    && expression.percentileCancer() < 0.05) {
                result.add(expression);
            }
        }

        return result;
    }

    @NotNull
    private static List<RnaFusion> selectNovelKnownFusions(@NotNull List<RnaFusion> rnaFusions, @NotNull List<LinxFusion> linxFusions,
            @NotNull KnownFusionCache knownFusionCache) {
        List<RnaFusion> result = Lists.newArrayList();

        for (RnaFusion rnaFusion : rnaFusions) {
            String geneUp = geneUp(rnaFusion);
            String geneDown = geneDown(rnaFusion);
            if (geneUp != null && geneDown != null) {
                if (knownFusionCache.hasKnownFusion(geneUp, geneDown) && !hasReportedLinxFusion(linxFusions, geneUp, geneDown)) {
                    result.add(rnaFusion);
                }
            }
        }

        return result;
    }

    @NotNull
    private static List<RnaFusion> selectNovelPromiscuousFusions(@NotNull List<RnaFusion> rnaFusions, @NotNull List<LinxFusion> linxFusions,
            @NotNull KnownFusionCache knownFusionCache) {
        List<RnaFusion> result = Lists.newArrayList();

        for (RnaFusion rnaFusion : rnaFusions) {
            boolean isTypeMatch = rnaFusion.svType().equals("BND") || rnaFusion.svType().equals("INV") || rnaFusion.svType().equals("INS");
            boolean hasSufficientDistance = Math.abs(rnaFusion.positionUp() - rnaFusion.positionDown()) > 1E6;

            String geneUp = geneUp(rnaFusion);
            String geneDown = geneDown(rnaFusion);
            if (geneUp != null && geneDown != null && (isTypeMatch || hasSufficientDistance)) {
                boolean isPromiscuous =
                        knownFusionCache.hasPromiscuousFiveGene(geneUp) || knownFusionCache.hasPromiscuousThreeGene(geneDown);
                boolean isKnown = knownFusionCache.hasKnownFusion(geneUp, geneDown);
                if (isPromiscuous && !isKnown && !hasReportedLinxFusion(linxFusions, geneUp, geneDown)) {
                    result.add(rnaFusion);
                }
            }
        }

        return result;
    }

    @VisibleForTesting
    @Nullable
    static String geneUp(@NotNull RnaFusion rnaFusion) {
        int split = rnaFusion.name().indexOf(RNA_FUSION_NAME_DELIMITER);
        return split > 0 ? rnaFusion.name().substring(0, split) : null;
    }

    @VisibleForTesting
    @Nullable
    static String geneDown(@NotNull RnaFusion rnaFusion) {
        int split = rnaFusion.name().indexOf(RNA_FUSION_NAME_DELIMITER);
        return split >= 0 && split < rnaFusion.name().length() - 1 ? rnaFusion.name().substring(split + 1) : null;
    }

    private static boolean hasReportedLinxFusion(@NotNull List<LinxFusion> linxFusions, @NotNull String geneUp, @NotNull String geneDown) {
        for (LinxFusion linxFusion : linxFusions) {
            if (linxFusion.reported() && linxFusion.geneStart().equals(geneUp) && linxFusion.geneEnd().equals(geneDown)) {
                return true;
            }
        }
        return false;
    }

    @NotNull
    private static List<NovelSpliceJunction> selectSkippedExons(@NotNull List<NovelSpliceJunction> junctions,
            @NotNull List<LinxFusion> linxFusions, @NotNull KnownFusionCache knownFusionCache) {
        List<NovelSpliceJunction> result = Lists.newArrayList();

        for (NovelSpliceJunction junction : junctions) {
            if (knownFusionCache.hasExonDelDup(junction.geneName())) {
                boolean isTypeMatch = junction.type().equals("SKIPPED_EXONS");
                boolean hasSufficientFragments = junction.fragmentCount() > 5;
                boolean hasLimitedCohortFreq = junction.cohortFrequency() < 30;
                boolean hasReportedLinxFusion = hasReportedLinxFusion(linxFusions, junction.geneName(), junction.geneName());
                if (isTypeMatch && hasSufficientFragments && hasLimitedCohortFreq && !hasReportedLinxFusion) {
                    result.add(junction);
                }
            }
        }

        return result;
    }

    @NotNull
    private static List<NovelSpliceJunction> selectNovelExonsIntrons(@NotNull List<NovelSpliceJunction> junctions,
            @NotNull List<DriverGene> driverGenes) {
        List<NovelSpliceJunction> result = Lists.newArrayList();

        Set<String> drivers = Sets.newHashSet();
        for (DriverGene driverGene : driverGenes) {
            drivers.add(driverGene.gene());
        }

        for (NovelSpliceJunction junction : junctions) {
            if (drivers.contains(junction.geneName())) {
                boolean isTypeMatch = junction.type().equals("NOVEL_INTRON") || junction.type().equals("NOVEL_EXON");
                boolean hasSufficientFragments = junction.fragmentCount() > 5;
                boolean hasLimitedCohortFreq = junction.cohortFrequency() < 10;
                if (isTypeMatch && hasSufficientFragments && hasLimitedCohortFreq) {
                    result.add(junction);
                }
            }
        }

        return result;
    }

    @NotNull
    private static Set<String> extractGenesOfType(@NotNull List<DriverGene> driverGenes, @NotNull DriverCategory categoryToInclude) {
        Set<String> filtered = Sets.newHashSet();
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.likelihoodType() == categoryToInclude) {
                filtered.add(driverGene.gene());
            }
        }
        return filtered;
    }
}
