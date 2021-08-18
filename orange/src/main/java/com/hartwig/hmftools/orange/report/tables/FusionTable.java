package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;

public final class FusionTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");

    private FusionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<LinxFusion> fusions) {
        if (fusions.isEmpty()) {
            return TableUtil.createEmptyTable(title, width);
        }

        Table table = TableUtil.createReportContentTable(width,
                new float[] { 1, 5 },
                new Cell[] { TableUtil.createHeaderCell("Fusion"), TableUtil.createHeaderCell("Details") });

        for (LinxFusion fusion : sort(fusions)) {
            table.addCell(TableUtil.createContentCell(fusion.name()));

            Table details = new Table(UnitValue.createPercentArray(new float[] { 1, 3 }));
            details.addCell(TableUtil.createKeyCell("5' End"));
            details.addCell(TableUtil.createValueCell(
                    fusion.geneStart() + " " + fusion.geneContextStart() + " (" + fusion.geneTranscriptStart() + ")"));
            details.addCell(TableUtil.createKeyCell("3' Start"));
            details.addCell(TableUtil.createValueCell(
                    fusion.geneEnd() + " " + fusion.geneContextEnd() + " (" + fusion.geneTranscriptEnd() + ")"));
            details.addCell(TableUtil.createKeyCell("Junction CN"));
            details.addCell(TableUtil.createValueCell(SINGLE_DIGIT.format(fusion.junctionCopyNumber())));
            details.addCell(TableUtil.createKeyCell("Phasing"));
            details.addCell(TableUtil.createValueCell(fusion.phased().display()));
            details.addCell(TableUtil.createKeyCell("Reported type (DL)"));
            details.addCell(TableUtil.createValueCell(fusion.reportedType() + " (" + fusion.likelihood().display() + ")"));
            details.addCell(TableUtil.createKeyCell("Chain links (terminated?)"));
            details.addCell(TableUtil.createValueCell(fusion.chainLinks() + (fusion.chainTerminated() ? " (Yes)" : " (No)")));
            details.addCell(TableUtil.createKeyCell("Domains kept"));
            details.addCell(TableUtil.createValueCell(!fusion.domainsKept().isEmpty() ? fusion.domainsKept() : "-"));
            details.addCell(TableUtil.createKeyCell("Domains lost"));
            details.addCell(TableUtil.createValueCell(!fusion.domainsLost().isEmpty() ? fusion.domainsLost() : "-"));
            table.addCell(TableUtil.createContentCell(details));
        }

        return TableUtil.createWrappingReportTable(table, title);
    }

    @NotNull
    private static List<LinxFusion> sort(@NotNull List<LinxFusion> fusions) {
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
}
