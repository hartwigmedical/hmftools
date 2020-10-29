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
import com.hartwig.hmftools.protect.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BachelorDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(BachelorDataLoader.class);

    @NotNull
    public static BachelorData load(@NotNull String bachelorTsv, @NotNull PurpleData purpleData, @NotNull LinxData linxData)
            throws IOException {

        List<ReportableGermlineVariant> germlineVariants = ReportableGermlineVariantFile.read(bachelorTsv);

        Set<String> somaticGenes = Sets.newHashSet();
        linxData.homozygousDisruptions().stream().map(ReportableHomozygousDisruption::gene).forEach(somaticGenes::add);
        linxData.geneDisruptions().stream().map(ReportableGeneDisruption::gene).forEach(somaticGenes::add);
        purpleData.somaticVariants().stream().map(ReportableVariant::gene).forEach(somaticGenes::add);
        purpleData.copyNumberAlterations()
                .stream()
                .filter(x -> !x.interpretation().equals(CopyNumberInterpretation.GAIN))
                .map(ReportableGainLoss::gene)
                .forEach(somaticGenes::add);

        List<ReportableVariant> reportableVariants = ReportableVariantFactory.reportableGermlineVariants(somaticGenes, germlineVariants);

        LOGGER.info("Loaded BACHELOR data from {}", new File(bachelorTsv).getParent());
        LOGGER.info(" Reportable germline variants: {}", reportableVariants.size());
        return ImmutableBachelorData.builder().addAllGermlineVariants(reportableVariants).build();
    }
}
