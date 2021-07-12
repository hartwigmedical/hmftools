package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.GeneCopyNumberTable;
import com.hartwig.hmftools.orange.report.tables.SomaticVariantTable;
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

        addSomaticVariants(document);
        addSomaticAmpDels(document);

        document.add(new Paragraph("TODO: Add Somatic Disruptions").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Somatic Fusions").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Somatic Viral Presence").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add LINX Visualisations").addStyle(ReportResources.tableContentStyle()));
    }

    private void addSomaticVariants(@NotNull Document document) {
        String driverVariantTableTitle = "Driver Variants (" + report.purple().reportableSomaticVariants().size() + ")";
        Table driverVariantTable = SomaticVariantTable.build(driverVariantTableTitle, report.purple().reportableSomaticVariants());
        DocumentUtil.addCheckedTable(document, driverVariantTableTitle, driverVariantTable);

        List<ReportableVariant> nonDriverVariants = extractInterestingNonDrivers(report.purple().unreportedSomaticVariants());
        String nonDriverVariantTableTitle = "Other potentially relevant coding variants (" + nonDriverVariants.size() + ")";
        Table nonDriverVariantTable = SomaticVariantTable.build(nonDriverVariantTableTitle, nonDriverVariants);
        DocumentUtil.addCheckedTable(document, nonDriverVariantTableTitle, nonDriverVariantTable);
    }

    private void addSomaticAmpDels(@NotNull Document document) {
        String driverAmpsDelsTitle = "Driver amps/dels (" + report.purple().reportableGainsLosses().size() + ")";
        Table driverAmpsDelsTable = GeneCopyNumberTable.build(driverAmpsDelsTitle, report.purple().reportableGainsLosses());
        DocumentUtil.addCheckedTable(document, driverAmpsDelsTitle, driverAmpsDelsTable);

        List<ReportableGainLoss> sortedGains = sortedGains(report.purple().unreportedGainsLosses());
        String sortedGainsTitle = "Non-driver amps (" + sortedGains.size() + ")";
        Table sortedGainsTable = GeneCopyNumberTable.build(sortedGainsTitle, sortedGains.subList(0, Math.min(10, sortedGains.size())));
        DocumentUtil.addCheckedTable(document, sortedGainsTitle, sortedGainsTable);

        List<ReportableGainLoss> lossesNoAllosomes = lossesNoAllosomes(report.purple().unreportedGainsLosses());
        String unreportedLossesTitle = "Non-driver dels on autosomes (" + lossesNoAllosomes.size() + ")";
        Table unreportedLossesTable =
                GeneCopyNumberTable.build(unreportedLossesTitle, lossesNoAllosomes.subList(0, Math.min(10, lossesNoAllosomes.size())));
        DocumentUtil.addCheckedTable(document, unreportedLossesTitle, unreportedLossesTable);
    }

    @NotNull
    private static List<ReportableGainLoss> sortedGains(@NotNull List<ReportableGainLoss> unreportedGainsLosses) {
        List<ReportableGainLoss> gains = Lists.newArrayList();
        for (ReportableGainLoss gainLoss : unreportedGainsLosses) {
            if (gainLoss.interpretation() == CopyNumberInterpretation.FULL_GAIN
                    || gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_GAIN) {
                gains.add(gainLoss);
            }
        }

        return gains.stream().sorted((o1, o2) -> (int) (o1.copies() - o2.copies())).collect(Collectors.toList());
    }

    @NotNull
    private static List<ReportableGainLoss> lossesNoAllosomes(@NotNull List<ReportableGainLoss> unreportedGainsLosses) {
        List<ReportableGainLoss> lossesNoAllosomes = Lists.newArrayList();
        for (ReportableGainLoss gainLoss : unreportedGainsLosses) {
            if (!HumanChromosome.fromString(gainLoss.chromosome()).isAllosome() && (
                    gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS
                            || gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS)) {
                lossesNoAllosomes.add(gainLoss);
            }
        }
        return lossesNoAllosomes;
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

    @NotNull
    private static ReportableVariant toReportable(@NotNull SomaticVariant variant) {
        return ReportableVariantFactory.fromVariant(variant, ReportableVariantSource.SOMATIC).driverLikelihood(Double.NaN).build();
    }
}
