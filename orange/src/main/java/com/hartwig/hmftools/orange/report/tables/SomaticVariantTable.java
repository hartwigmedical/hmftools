package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SomaticVariantTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");
    private static final DecimalFormat PERCENTAGE_FORMAT = ReportResources.decimalFormat("#'%'");

    private SomaticVariantTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableVariant> variants) {
        if (variants.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 3, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Variant"), Cells.createHeader("VCN"), Cells.createHeader("CN"), Cells.createHeader("MACN"),
                        Cells.createHeader("RNA VAF"), Cells.createHeader("Biallelic"), Cells.createHeader("Hotspot"),
                        Cells.createHeader("DL"), Cells.createHeader("CL"), Cells.createHeader("Phase ID") });

        for (ReportableVariant variant : VariantUtil.sort(variants)) {
            table.addCell(Cells.createContent(VariantUtil.variantField(variant)));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(variant.alleleCopyNumber())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(variant.totalCopyNumber())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(variant.minorAlleleCopyNumber())));
            table.addCell(Cells.createContent(ReportResources.NOT_AVAILABLE));
            table.addCell(Cells.createContent(variant.biallelic() ? "Yes" : "No"));
            table.addCell(Cells.createContent(VariantUtil.hotspotField(variant)));
            table.addCell(Cells.createContent(driverLikelihoodField(variant)));
            table.addCell(Cells.createContent(PERCENTAGE_FORMAT.format(variant.clonalLikelihood() * 100)));
            table.addCell(Cells.createContent(phaseSetField(variant)));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static String driverLikelihoodField(@NotNull ReportableVariant variant) {
        if (Double.isNaN(variant.driverLikelihood())) {
            return Strings.EMPTY;
        } else {
            return PERCENTAGE_FORMAT.format(variant.driverLikelihood() * 100);
        }
    }

    @NotNull
    private static String phaseSetField(@NotNull ReportableVariant variant) {
        return variant.localPhaseSet() != null ? String.valueOf(variant.localPhaseSet()) : Strings.EMPTY;
    }
}
