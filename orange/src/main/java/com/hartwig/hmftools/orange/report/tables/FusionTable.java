package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.CellUtil;
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
            return TableUtil.createEmpty(title, width);
        }

        Table table = TableUtil.createContent(width,
                new float[] { 1, 5 },
                new Cell[] { CellUtil.createHeader("Fusion"), CellUtil.createHeader("Details") });

        for (LinxFusion fusion : sort(fusions)) {
            table.addCell(CellUtil.createContent(fusion.name()));

            Table details = new Table(UnitValue.createPercentArray(new float[] { 1, 3 }));
            details.addCell(CellUtil.createKey("5' End"));
            details.addCell(CellUtil.createValue(fiveEndString(fusion)));
            details.addCell(CellUtil.createKey("3' Start"));
            details.addCell(CellUtil.createValue(threeStartString(fusion)));
            details.addCell(CellUtil.createKey("Junction CN"));
            details.addCell(CellUtil.createValue(SINGLE_DIGIT.format(fusion.junctionCopyNumber())));
            details.addCell(CellUtil.createKey("RNA expression"));
            details.addCell(CellUtil.createValue(ReportResources.NOT_AVAILABLE));
            details.addCell(CellUtil.createKey("Phasing"));
            details.addCell(CellUtil.createValue(fusion.phased().displayStr()));
            details.addCell(CellUtil.createKey("Reported type (DL)"));
            details.addCell(CellUtil.createValue(fusion.reportedType() + " (" + fusion.likelihood().displayStr() + ")"));
            details.addCell(CellUtil.createKey("Chain links (terminated?)"));
            details.addCell(CellUtil.createValue(fusion.chainLinks() + (fusion.chainTerminated() ? " (Yes)" : " (No)")));
            details.addCell(CellUtil.createKey("Domains kept"));
            details.addCell(CellUtil.createValue(!fusion.domainsKept().isEmpty() ? fusion.domainsKept() : "-"));
            details.addCell(CellUtil.createKey("Domains lost"));
            details.addCell(CellUtil.createValue(!fusion.domainsLost().isEmpty() ? fusion.domainsLost() : "-"));
            // Need to keep this details table to avoid page-wrapping that cuts through the middle of a single fusion
            table.addCell(CellUtil.createContent(details).setKeepTogether(true));
        }

        return TableUtil.createWrapping(table, title);
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
