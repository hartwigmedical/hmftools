package com.hartwig.hmftools.protect.bachelor;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariantFile;
import com.hartwig.hmftools.protect.ProtectConfig;
import com.hartwig.hmftools.protect.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantFactory;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingFile;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class BachelorDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(BachelorDataLoader.class);

    private BachelorDataLoader() {
    }

    @NotNull
    public static BachelorData load(@NotNull ProtectConfig config, @NotNull PurpleData purpleData, @NotNull LinxData linxData)
            throws IOException {
        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromTsv(config.germlineReportingTsv());

        return load(config.bachelorTsv(), purpleData, linxData, germlineReportingModel);
    }

    @NotNull
    public static BachelorData load(@NotNull String bachelorTsv, @NotNull PurpleData purpleData, @NotNull LinxData linxData,
            @NotNull GermlineReportingModel germlineReportingModel) throws IOException {
        LOGGER.info("Loading BACHELOR data from {}", new File(bachelorTsv).getParent());
        List<ReportableGermlineVariant> germlineVariants = ReportableGermlineVariantFile.read(bachelorTsv);

        Set<String> genesWithInactivationEvent = Sets.newHashSet();
        linxData.homozygousDisruptions().stream().map(ReportableHomozygousDisruption::gene).forEach(genesWithInactivationEvent::add);
        linxData.geneDisruptions().stream().map(ReportableGeneDisruption::gene).forEach(genesWithInactivationEvent::add);
        purpleData.somaticVariants().stream().map(ReportableVariant::gene).forEach(genesWithInactivationEvent::add);
        purpleData.copyNumberAlterations()
                .stream()
                .filter(x -> !x.interpretation().equals(CopyNumberInterpretation.GAIN))
                .map(ReportableGainLoss::gene)
                .forEach(genesWithInactivationEvent::add);

        List<ReportableVariant> reportableVariants = ReportableVariantFactory.reportableGermlineVariants(germlineVariants,
                genesWithInactivationEvent,
                germlineReportingModel);

        LOGGER.info(" Loaded {} reportable germline variants from {}", reportableVariants.size(), bachelorTsv);
        return ImmutableBachelorData.builder().addAllGermlineVariants(reportableVariants).build();
    }
}
