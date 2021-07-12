package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.sv.linx.FusionPhasedType;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.FusionTable;
import com.hartwig.hmftools.orange.report.tables.GeneCopyNumberTable;
import com.hartwig.hmftools.orange.report.tables.GeneDisruptionTable;
import com.hartwig.hmftools.orange.report.tables.SomaticVariantTable;
import com.hartwig.hmftools.orange.report.tables.ViralPresenceTable;
import com.hartwig.hmftools.orange.report.util.ImageUtil;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SomaticFindingsChapter implements ReportChapter {

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
        addFusions(document);
        addViralPresence(document);
        addSomaticGeneDisruptions(document);
        addLinxDriverPlots(document);
    }

    private void addSomaticVariants(@NotNull Document document) {
        String driverVariantTableTitle = "Driver variants (" + report.purple().reportableSomaticVariants().size() + ")";
        Table driverVariantTable = SomaticVariantTable.build(driverVariantTableTitle, report.purple().reportableSomaticVariants());
        document.add(driverVariantTable);

        List<ReportableVariant> nonDriverVariants = extractInterestingNonDrivers(report.purple().unreportedSomaticVariants());
        String nonDriverVariantTableTitle = "Other potentially relevant coding variants (" + nonDriverVariants.size() + ")";
        Table nonDriverVariantTable = SomaticVariantTable.build(nonDriverVariantTableTitle, nonDriverVariants);
        document.add(nonDriverVariantTable);
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

    private void addSomaticAmpDels(@NotNull Document document) {
        String driverAmpsDelsTitle = "Driver amps/dels (" + report.purple().reportableGainsLosses().size() + ")";
        Table driverAmpsDelsTable = GeneCopyNumberTable.build(driverAmpsDelsTitle, report.purple().reportableGainsLosses());
        document.add(driverAmpsDelsTable);

        List<ReportableGainLoss> sortedGains = sortedGains(report.purple().unreportedGainsLosses());
        String sortedGainsTitle = "Other amps (" + sortedGains.size() + ")";
        Table sortedGainsTable = GeneCopyNumberTable.build(sortedGainsTitle, sortedGains.subList(0, Math.min(10, sortedGains.size())));
        document.add(sortedGainsTable);

        List<ReportableGainLoss> lossesNoAllosomes = lossesNoAllosomes(report.purple().unreportedGainsLosses());
        String unreportedLossesTitle = "Other dels on autosomes (" + lossesNoAllosomes.size() + ")";
        Table unreportedLossesTable =
                GeneCopyNumberTable.build(unreportedLossesTitle, lossesNoAllosomes.subList(0, Math.min(10, lossesNoAllosomes.size())));
        document.add(unreportedLossesTable);
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

    private void addFusions(@NotNull Document document) {
        String driverFusionTitle = "Driver fusions (" + report.linx().reportableFusions().size() + ")";
        Table driverFusionTable = FusionTable.build(driverFusionTitle, report.linx().reportableFusions());
        document.add(driverFusionTable);

        List<LinxFusion> inframes = inframe(report.linx().unreportedFusions());
        String inframeNonDriverTitle = "Other inframe fusions (" + inframes.size() + ")";
        Table inframeNonDriverTable = FusionTable.build(inframeNonDriverTitle, inframes.subList(0, Math.min(10, inframes.size())));
        document.add(inframeNonDriverTable);
    }

    @NotNull
    private static List<LinxFusion> inframe(@NotNull List<LinxFusion> fusions) {
        List<LinxFusion> filtered = Lists.newArrayList();
        for (LinxFusion fusion : fusions) {
            if (fusion.phased() == FusionPhasedType.INFRAME) {
                filtered.add(fusion);
            }
        }
        return filtered;
    }

    private void addViralPresence(@NotNull Document document) {
        String driverVirusTitle = "Driver viruses (" + report.virusInterpreter().reportableViruses().size() + ")";
        Table driverVirusTable = ViralPresenceTable.build(driverVirusTitle, report.virusInterpreter().reportableViruses());
        document.add(driverVirusTable);

        String otherVirusTitle = "Other viral presence (" + report.virusInterpreter().unreportedViruses().size() + ")";
        Table otherVirusTable = ViralPresenceTable.build(otherVirusTitle, report.virusInterpreter().unreportedViruses());
        document.add(otherVirusTable);
    }

    private void addSomaticGeneDisruptions(@NotNull Document document) {
        String geneDisruptionTitle = "Gene disruptions (" + report.linx().geneDisruptions().size() + ")";
        Table geneDisruptionTable = GeneDisruptionTable.build(geneDisruptionTitle, report.linx().geneDisruptions());
        document.add(geneDisruptionTable);
    }

    private void addLinxDriverPlots(@NotNull Document document) {
        document.add(new Paragraph("Linx driver plots").addStyle(ReportResources.tableTitleStyle()));
        for (String plot : report.plots().linxDriverPlots()) {
            Image image = ImageUtil.build(plot);
            image.setMaxHeight(300);
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            document.add(image);
        }
    }
}
