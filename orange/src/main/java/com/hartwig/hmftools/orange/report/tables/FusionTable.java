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

import org.jetbrains.annotations.NotNull;

public final class FusionTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");

    private FusionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, @NotNull List<LinxFusion> fusions) {
        if (fusions.isEmpty()) {
            return TableUtil.createEmptyTable(title);
        }

        Table table = TableUtil.createReportContentTable(new float[] { 80, 80, 80, 40, 40, 40, 65, 40 },
                new Cell[] { TableUtil.createHeaderCell("Fusion"), TableUtil.createHeaderCell("5' Transcript"),
                        TableUtil.createHeaderCell("3' Transcript"), TableUtil.createHeaderCell("5' End"),
                        TableUtil.createHeaderCell("3' Start"), TableUtil.createHeaderCell("Copies"),
                        TableUtil.createHeaderCell("Phasing"), TableUtil.createHeaderCell("Driver") });

        for (LinxFusion fusion : sort(fusions)) {
            table.addCell(TableUtil.createContentCell(fusion.name()));
            table.addCell(TableUtil.createContentCell(fusion.geneTranscriptStart()));
            table.addCell(TableUtil.createContentCell(fusion.geneTranscriptEnd()));
            table.addCell(TableUtil.createContentCell(fusion.geneContextStart()));
            table.addCell(TableUtil.createContentCell(fusion.geneContextEnd()));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(fusion.junctionCopyNumber())));
            table.addCell(TableUtil.createContentCell(fusion.phased().display()));
            table.addCell(TableUtil.createContentCell(fusion.likelihood().display()));
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
