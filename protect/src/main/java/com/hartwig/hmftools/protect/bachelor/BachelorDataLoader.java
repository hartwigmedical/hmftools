package com.hartwig.hmftools.protect.bachelor;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariantFile;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantFactory;
import com.hartwig.hmftools.protect.variants.germline.DriverGermlineVariant;
import com.hartwig.hmftools.protect.variants.germline.FilterGermlineVariants;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BachelorDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(BachelorDataLoader.class);

    public static BachelorData load(@NotNull PurpleData purpleData, @NotNull String bachelorTsv, @NotNull ChordStatus chordStatus)
            throws IOException {
        List<ReportableVariant> result = Lists.newArrayList();

        List<ReportableGermlineVariant> variants = ReportableGermlineVariantFile.read(bachelorTsv);

        // Every gene is allowed at this stage.
        Set<String> allowedGenes = variants.stream().map(ReportableGermlineVariant::gene).collect(Collectors.toSet());

        List<DriverGermlineVariant> driverGermlineVariants = FilterGermlineVariants.filterGermlineVariantsForReporting(variants,
                allowedGenes,
                purpleData.geneCopyNumbers(),
                purpleData.somaticVariants().stream().map(ReportableVariant::gene).collect(Collectors.toSet()),
                chordStatus);

        result.addAll(ReportableVariantFactory.reportableGermlineVariants(driverGermlineVariants));

        LOGGER.info("Loaded BACHELOR data from {}", new File(bachelorTsv).getParent());
        LOGGER.info(" Reportable germline variants: {}", result.size());
        return ImmutableBachelorData.builder().addAllGermlineVariants(result).build();

    }

}
