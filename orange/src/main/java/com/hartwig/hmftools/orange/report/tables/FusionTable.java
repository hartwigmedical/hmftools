package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.isofox.IsofoxData;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.IBlockElement;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FusionTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");
    private static final DecimalFormat PERCENTAGE = ReportResources.decimalFormat("#'%'");

    private FusionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<LinxFusion> fusions, @Nullable IsofoxData isofox) {
        if (fusions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 5 },
                new Cell[] { Cells.createHeader("Fusion"), Cells.createHeader("Details") });

        for (LinxFusion fusion : sortLinxFusions(fusions)) {
            table.addCell(Cells.createContent(fusion.name()));

            Table details = new Table(UnitValue.createPercentArray(new float[] { 1, 3 }));
            details.addCell(Cells.createKey("5' End"));
            details.addCell(Cells.createValue(fiveEndString(fusion)));
            details.addCell(Cells.createKey("3' Start"));
            details.addCell(Cells.createValue(threeStartString(fusion)));
            details.addCell(Cells.createKey("Junction CN"));
            details.addCell(Cells.createValue(SINGLE_DIGIT.format(fusion.junctionCopyNumber())));
            details.addCell(Cells.createKey("RNA support"));
            details.addCell(Cells.createValue(rnaFragmentSupportTable(isofox, fusion)).setKeepTogether(true));
            details.addCell(Cells.createKey("Phasing"));
            details.addCell(Cells.createValue(fusion.phased().displayStr()));
            details.addCell(Cells.createKey("Reported type (DL)"));
            details.addCell(Cells.createValue(fusion.reportedType() + " (" + fusion.likelihood().displayStr() + ")"));
            details.addCell(Cells.createKey("Chain links (terminated?)"));
            details.addCell(Cells.createValue(fusion.chainLinks() + (fusion.chainTerminated() ? " (Yes)" : " (No)")));
            details.addCell(Cells.createKey("Domains kept"));
            details.addCell(Cells.createValue(!fusion.domainsKept().isEmpty() ? fusion.domainsKept() : "-"));
            details.addCell(Cells.createKey("Domains lost"));
            details.addCell(Cells.createValue(!fusion.domainsLost().isEmpty() ? fusion.domainsLost() : "-"));
            // Need to keep this details table to avoid page-wrapping that cuts through the middle of a single fusion
            table.addCell(Cells.createContent(details).setKeepTogether(true));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static String fiveEndString(@NotNull LinxFusion fusion) {
        return fusion.geneStart() + " " + fusion.geneContextStart() + " (" + fusion.geneTranscriptStart() + ")";
    }

    @NotNull
    private static String threeStartString(@NotNull LinxFusion fusion) {
        return fusion.geneEnd() + " " + fusion.geneContextEnd() + " (" + fusion.geneTranscriptEnd() + ")";
    }

    @NotNull
    private static IBlockElement rnaFragmentSupportTable(@Nullable IsofoxData isofox, @NotNull LinxFusion fusion) {
        if (isofox == null) {
            return new Paragraph(ReportResources.NOT_AVAILABLE);
        }

        if (fusion.reportedType().equals(KnownFusionType.IG_KNOWN_PAIR.toString()) || fusion.reportedType()
                .equals(KnownFusionType.IG_PROMISCUOUS.toString())) {
            return supportFromExpressionOfGeneEnd(isofox, fusion);
        } else if (fusion.reportedType().equals(KnownFusionType.EXON_DEL_DUP.toString())) {
            return supportFromSpliceJunctions(isofox, fusion);
        } else {
            return supportFromRnaFusions(isofox, fusion);
        }
    }

    @NotNull
    private static IBlockElement supportFromExpressionOfGeneEnd(@NotNull IsofoxData isofox, @NotNull LinxFusion fusion) {
        GeneExpression threeExpression = null;
        for (GeneExpression geneExpression : isofox.geneExpressions()) {
            if (geneExpression.geneName().equals(fusion.geneEnd())) {
                threeExpression = geneExpression;
                break;
            }
        }

        if (threeExpression == null) {
            return new Paragraph("None");
        }

        String tpmString = "TPM " + SINGLE_DIGIT.format(threeExpression.tpm());
        String fcTypeString = "FC " + SINGLE_DIGIT.format(threeExpression.tpm() / threeExpression.medianTpmCancer());
        String typeString = " Type percentile " + PERCENTAGE.format(threeExpression.percentileCancer() * 100) + " (" + fcTypeString + ")";
        String fcDbString = "FC " + SINGLE_DIGIT.format(threeExpression.tpm() / threeExpression.medianTpmCohort());
        String dbString = " DB percentile " + PERCENTAGE.format(threeExpression.percentileCohort() * 100) + " (" + fcDbString + ")";

        return new Paragraph(threeExpression.geneName() + " " + tpmString + ", " + typeString + ", " + dbString);
    }

    @NotNull
    private static IBlockElement supportFromSpliceJunctions(@NotNull IsofoxData isofox, @NotNull LinxFusion fusion) {
        List<NovelSpliceJunction> matches = Lists.newArrayList();
        for (NovelSpliceJunction junction : isofox.novelSpliceJunctions()) {
            if (junction.geneName().equals(fusion.geneStart()) && junction.geneName().equals(fusion.geneEnd())) {
                matches.add(junction);
            }
        }

        if (matches.isEmpty()) {
            return new Paragraph("None");
        }

        Table fragmentSupportTable = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        for (NovelSpliceJunction junction : max5(sortNovelSpliceJunctions(matches))) {
            fragmentSupportTable.addCell(Cells.createKey(
                    junction.chromosome() + ":" + junction.junctionStart() + "-" + junction.junctionEnd()));

            String fragments = junction.fragmentCount() + " fragments";
            String depth = junction.depthStart() + " / " + junction.depthEnd() + " depth";
            fragmentSupportTable.addCell(Cells.createValue(junction.type() + " (" + fragments + ", " + depth + ")"));
        }

        return fragmentSupportTable;
    }

    @NotNull
    private static List<NovelSpliceJunction> sortNovelSpliceJunctions(@NotNull List<NovelSpliceJunction> novelSpliceJunctions) {
        return novelSpliceJunctions.stream()
                .sorted((junction1, junction2) -> Integer.compare(junction2.fragmentCount(), junction1.fragmentCount()))
                .collect(Collectors.toList());
    }

    @NotNull
    private static IBlockElement supportFromRnaFusions(@NotNull IsofoxData isofox, @NotNull LinxFusion fusion) {
        List<RnaFusion> matches = Lists.newArrayList();
        for (RnaFusion rnaFusion : isofox.fusions()) {
            if (rnaFusion.name().equals(fusion.name())) {
                matches.add(rnaFusion);
            }
        }

        if (matches.isEmpty()) {
            return new Paragraph("None");
        }

        Table fragmentSupportTable = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        for (RnaFusion rnaFusion : max5(sortRnaFusions(matches))) {
            String up = rnaFusion.chromosomeUp() + ":" + rnaFusion.positionUp();
            String down = rnaFusion.chromosomeDown() + ":" + rnaFusion.positionDown();
            fragmentSupportTable.addCell(Cells.createKey(up + "-" + down));

            String split = rnaFusion.splitFragments() + " split";
            String realigned = rnaFusion.realignedFrags() + " realig.";
            String discord = rnaFusion.discordantFrags() + " discord.";
            String fragments = split + " / " + realigned + " / " + discord + " fragments";
            String depth = rnaFusion.depthUp() + " / " + rnaFusion.depthDown() + " depth";
            fragmentSupportTable.addCell(Cells.createValue(fragments + ", " + depth));
        }

        return fragmentSupportTable;
    }

    @NotNull
    private static List<RnaFusion> sortRnaFusions(@NotNull List<RnaFusion> rnaFusions) {
        return rnaFusions.stream().sorted((fusion1, fusion2) -> {
            int sumFragments1 = fusion1.splitFragments() + fusion1.realignedFrags() + fusion1.discordantFrags();
            int sumFragments2 = fusion2.splitFragments() + fusion2.realignedFrags() + fusion2.discordantFrags();
            return Integer.compare(sumFragments2, sumFragments1);
        }).collect(Collectors.toList());
    }

    @NotNull
    private static List<LinxFusion> sortLinxFusions(@NotNull List<LinxFusion> fusions) {
        return fusions.stream().sorted((fusion1, fusion2) -> {
            if (fusion1.likelihood() == fusion2.likelihood()) {
                if (fusion1.geneStart().equals(fusion2.geneStart())) {
                    return fusion1.geneEnd().compareTo(fusion2.geneEnd());
                } else {
                    return fusion1.geneStart().compareTo(fusion2.geneStart());
                }
            } else {
                return fusion1.likelihood() == FusionLikelihoodType.HIGH ? -1 : 1;
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static <T> List<T> max5(@NotNull List<T> elements) {
        return elements.subList(0, Math.min(5, elements.size()));
    }
}
