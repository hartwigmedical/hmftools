package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.SomaticDriverTable;
import com.hartwig.hmftools.orange.report.util.DocumentUtil;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SomaticFindingsChapter implements ReportChapter {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("0.0");
    private static final DecimalFormat PERCENTAGE_FORMAT = ReportResources.decimalFormat("#'%'");

    @NotNull
    private final OrangeReport report;

    public SomaticFindingsChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Somatic Findings";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Somatic Findings").addStyle(ReportResources.chapterTitleStyle()));

        Table driverVariantTable = SomaticDriverTable.build("Driver Variants (" + report.purple().reportableSomaticVariants().size() + ")",
                report.purple().reportableSomaticVariants());
        DocumentUtil.addCheckedTable(document, driverVariantTable, "No driver variants found");

        List<ReportableVariant> nonDriverVariants = extractInterestingNonDrivers(report.purple().unreportedSomaticVariants());
        Table nonDriverVariantTable =
                SomaticDriverTable.build("Other coding variants (" + report.purple().unreportedSomaticVariants().size() + ")",
                        nonDriverVariants);
        DocumentUtil.addCheckedTable(document, nonDriverVariantTable, "No interesting other coding variants found!");

        document.add(new Paragraph("TODO: Add Somatic AMPs/DELs").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Somatic Disruptions").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Somatic Fusions").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Somatic Viral Presence").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add LINX Visualisations").addStyle(ReportResources.tableContentStyle()));
    }

    @NotNull
    private static List<ReportableVariant> extractInterestingNonDrivers(@NotNull List<SomaticVariant> variants) {
        List<ReportableVariant> filtered = Lists.newArrayList();
        for (SomaticVariant variant : variants) {
            if (!variant.reported()) {
                if (variant.isHotspot() || hasReportedVariantWithPhase(variants, variant.localPhaseSet())) {
                    filtered.add(toReportable(variant));
                }
            }
        }
        return filtered;
    }

    private static boolean hasReportedVariantWithPhase(@NotNull List<SomaticVariant> variants, @Nullable Integer targetPhaseSet) {
        if (targetPhaseSet == null) {
            return false;
        }

        for (SomaticVariant variant : variants) {
            if (variant.reported() && variant.localPhaseSet() != null && variant.localPhaseSet().equals(targetPhaseSet)) {
                return true;
            }
        }

        return false;
    }

    private static ReportableVariant toReportable(@NotNull SomaticVariant variant) {
        return ReportableVariantFactory.fromVariant(variant, ReportableVariantSource.SOMATIC).driverLikelihood(Double.NaN).build();
    }
}
