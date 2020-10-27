package com.hartwig.hmftools.protect.bachelor;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariantFile;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantFactory;
import com.hartwig.hmftools.protect.variants.germline.DriverGermlineVariant;
import com.hartwig.hmftools.protect.variants.germline.FilterGermlineVariants;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BachelorDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(BachelorDataLoader.class);

    @NotNull
    public static BachelorData load(@NotNull String bachelorTsv, @NotNull PurpleData purpleData, @NotNull LinxData linxData)
            throws IOException {
        List<ReportableVariant> result = Lists.newArrayList();

        List<ReportableGermlineVariant> germlineVariants = ReportableGermlineVariantFile.read(bachelorTsv);

        // Every gene is allowed at this stage.
        Set<String> allowedGenes = germlineVariants.stream().map(ReportableGermlineVariant::gene).collect(Collectors.toSet());

        Set<String> genesWithSomaticInactivationEvent = Sets.newHashSet();
        genesWithSomaticInactivationEvent.addAll(purpleData.somaticVariants()
                .stream()
                .map(ReportableVariant::gene)
                .collect(Collectors.toSet()));

        genesWithSomaticInactivationEvent.addAll(linxData.geneDisruptions()
                .stream()
                .map(ReportableGeneDisruption::gene)
                .collect(Collectors.toSet()));
        // TODO: Complete list of genes with somatic variant (also hom disruptions and losses).

        List<DriverGermlineVariant> driverGermlineVariants = FilterGermlineVariants.filterGermlineVariantsForReporting(allowedGenes,
                germlineVariants,
                genesWithSomaticInactivationEvent);

        result.addAll(ReportableVariantFactory.reportableGermlineVariants(driverGermlineVariants));

        LOGGER.info("Loaded BACHELOR data from {}", new File(bachelorTsv).getParent());
        LOGGER.info(" Reportable germline variants: {}", result.size());
        return ImmutableBachelorData.builder().addAllGermlineVariants(result).build();
    }
}
