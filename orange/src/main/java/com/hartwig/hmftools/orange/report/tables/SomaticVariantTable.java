package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.CellUtil;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.hartwig.hmftools.orange.report.util.VariantUtil;
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
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableVariant> driverVariants) {
        if (driverVariants.isEmpty()) {
            return TableUtil.createEmpty(title, width);
        }

        Table table = TableUtil.createContent(width,
                new float[] { 3, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { CellUtil.createHeader("Variant"), CellUtil.createHeader("VCN"), CellUtil.createHeader("CN"),
                        CellUtil.createHeader("MACN"), CellUtil.createHeader("RNA VAF"), CellUtil.createHeader("Biallelic"),
                        CellUtil.createHeader("Hotspot"), CellUtil.createHeader("DL"), CellUtil.createHeader("CL"),
                        CellUtil.createHeader("Phase ID") });

        for (ReportableVariant variant : VariantUtil.sort(driverVariants)) {
            table.addCell(CellUtil.createContent(VariantUtil.variantField(variant)));
            table.addCell(CellUtil.createContent(SINGLE_DIGIT.format(variant.alleleCopyNumber())));
            table.addCell(CellUtil.createContent(SINGLE_DIGIT.format(variant.totalCopyNumber())));
            table.addCell(CellUtil.createContent(SINGLE_DIGIT.format(variant.minorAlleleCopyNumber())));
            table.addCell(CellUtil.createContent(ReportResources.NOT_AVAILABLE));
            table.addCell(CellUtil.createContent(variant.biallelic() ? "Yes" : "No"));
            table.addCell(CellUtil.createContent(VariantUtil.hotspotField(variant)));
            table.addCell(CellUtil.createContent(driverLikelihoodField(variant)));
            table.addCell(CellUtil.createContent(PERCENTAGE_FORMAT.format(variant.clonalLikelihood() * 100)));
            table.addCell(CellUtil.createContent(phaseSetField(variant)));
        }

        return TableUtil.createWrapping(table, title);
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
