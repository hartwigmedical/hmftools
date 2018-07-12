package com.hartwig.hmftools.healthchecker.runners;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.NoResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class StrelkaChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(StrelkaChecker.class);

    private static final String STRELKA_OUTPUT_EXTENSION = "_post_processed.vcf.gz";

    public StrelkaChecker() {
    }

    @NotNull
    @Override
    public BaseResult run(@NotNull final RunContext runContext) {
        if (!runContext.isSomaticRun()) {
            return new NoResult(CheckType.STRELKA);
        }
        final List<SomaticVariant> variants;
        try {
            final String vcfPath = vcfFilePath(runContext.runDirectory(), STRELKA_OUTPUT_EXTENSION);
            variants = SomaticVariantFactory.unfilteredInstance().fromVCFFile(runContext.tumorSample(), vcfPath);
        } catch (IOException exc) {
            LOGGER.warn("Could not load strelka variants.");
            return new NoResult(CheckType.STRELKA);
        }
        final List<SomaticVariant> snps = extractVariantsWithType(variants, VariantType.SNP);
        final List<SomaticVariant> mnps = extractVariantsWithType(variants, VariantType.MNP);
        final List<SomaticVariant> indels = extractVariantsWithType(variants, VariantType.INDEL);

        final List<SomaticVariant> snpsWithDBSNPAndNotCOSMIC = extractDBSNPNotCOSMIC(snps);
        final List<SomaticVariant> mnpWithDBSNPAndNotCOSMIC = extractDBSNPNotCOSMIC(mnps);
        final List<SomaticVariant> indelsWithDBSNPAndNotCOSMIC = extractDBSNPNotCOSMIC(indels);

        final List<HealthCheck> checks = Lists.newArrayList();
        checks.add(new HealthCheck(runContext.tumorSample(), StrelkaCheck.SOMATIC_SNP_COUNT.toString(), String.valueOf(snps.size())));
        checks.add(new HealthCheck(runContext.tumorSample(), StrelkaCheck.SOMATIC_INDEL_COUNT.toString(), String.valueOf(indels.size())));
        checks.add(new HealthCheck(runContext.tumorSample(), StrelkaCheck.SOMATIC_SNP_DBSNP_COUNT.toString(),
                String.valueOf(snpsWithDBSNPAndNotCOSMIC.size())));
        checks.add(new HealthCheck(runContext.tumorSample(), StrelkaCheck.SOMATIC_INDEL_DBSNP_COUNT.toString(),
                String.valueOf(indelsWithDBSNPAndNotCOSMIC.size())));
        checks.add(new HealthCheck(runContext.tumorSample(), StrelkaCheck.SOMATIC_MNP_COUNT.toString(), String.valueOf(mnps.size())));
        checks.add(new HealthCheck(runContext.tumorSample(), StrelkaCheck.SOMATIC_MNP_DBSNP_COUNT.toString(),
                String.valueOf(mnpWithDBSNPAndNotCOSMIC.size())));

        return toMultiValueResult(checks);
    }

    @NotNull
    private static String vcfFilePath(@NotNull String basePath, @NotNull String fileExtension) throws FileNotFoundException {
        return PathExtensionFinder.build().findPath(basePath, fileExtension).toString();
    }

    @NotNull
    private static BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);
        return new MultiValueResult(CheckType.STRELKA, checks);
    }

    @NotNull
    private static List<SomaticVariant> extractVariantsWithType(@NotNull List<SomaticVariant> variants, @NotNull VariantType type) {
        return filter(variants, variant -> variant.type().equals(type));
    }

    @NotNull
    private static List<SomaticVariant> extractDBSNPNotCOSMIC(@NotNull List<SomaticVariant> variants) {
        return filter(variants, variant -> !variant.isCOSMIC() && variant.isDBSNP());
    }

    @NotNull
    private static List<SomaticVariant> filter(@NotNull List<SomaticVariant> variants, @NotNull Predicate<SomaticVariant> filter) {
        return variants.stream().filter(filter).collect(Collectors.toList());
    }
}
