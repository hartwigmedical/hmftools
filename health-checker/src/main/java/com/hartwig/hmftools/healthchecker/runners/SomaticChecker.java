package com.hartwig.hmftools.healthchecker.runners;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantConstants;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.predicate.VariantPredicates;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.SomaticCheck;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.SOMATIC)
public class SomaticChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(SomaticChecker.class);

    private static final String SOMATICS_EXTENSION = "_melted.vcf";
    private static final List<Integer> CALLERS_COUNT = Arrays.asList(1, 2, 3, 4);

    private static final double AF_SD_DISTANCE = 0.16;

    public SomaticChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.SOMATIC;
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HartwigException {
        final VCFSomaticFile variantFile = VCFFileLoader.loadSomaticVCF(runContext.runDirectory(), SOMATICS_EXTENSION);

        if (!variantFile.sample().equals(runContext.tumorSample())) {
            LOGGER.warn("Sample name in VCF (" + variantFile.sample() + ") does not match with name ("
                    + runContext.tumorSample() + ") from run context!");
        }
        final List<SomaticVariant> variants = passOnly(variantFile.variants());

        final List<HealthCheck> checks = Lists.newArrayList();
        checks.addAll(getTypeChecks(variants, runContext.tumorSample(), VariantType.SNP));
        checks.addAll(getTypeChecks(variants, runContext.tumorSample(), VariantType.INDEL));
        checks.addAll(getAFChecks(variants, runContext.tumorSample()));

        return toMultiValueResult(checks);
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        final List<HealthCheck> checks = Lists.newArrayList();
        for (final VariantType type : VariantType.values()) {
            checks.add(new HealthCheck(runContext.tumorSample(), SomaticCheck.COUNT_TOTAL.checkName(type.name()),
                    HealthCheckConstants.ERROR_VALUE));
            checks.add(new HealthCheck(runContext.tumorSample(), SomaticCheck.DBSNP_COUNT.checkName(type.name()),
                    HealthCheckConstants.ERROR_VALUE));
            checks.addAll(SomaticVariantConstants.ALL_CALLERS.stream().map(
                    caller -> new HealthCheck(runContext.tumorSample(),
                            SomaticCheck.COUNT_PER_CALLER.checkName(type.name(), caller),
                            HealthCheckConstants.ERROR_VALUE)).collect(Collectors.toList()));
            checks.addAll(SomaticVariantConstants.ALL_CALLERS.stream().map(
                    caller -> new HealthCheck(runContext.tumorSample(),
                            SomaticCheck.PRECISION_CHECK.checkName(type.name(), caller),
                            HealthCheckConstants.ERROR_VALUE)).collect(Collectors.toList()));
            checks.addAll(SomaticVariantConstants.ALL_CALLERS.stream().map(
                    caller -> new HealthCheck(runContext.tumorSample(),
                            SomaticCheck.SENSITIVITY_CHECK.checkName(type.name(), caller),
                            HealthCheckConstants.ERROR_VALUE)).collect(Collectors.toList()));
            checks.addAll(CALLERS_COUNT.stream().map(count -> new HealthCheck(runContext.tumorSample(),
                    SomaticCheck.PROPORTION_CHECK.checkName(type.name(), String.valueOf(count)),
                    HealthCheckConstants.ERROR_VALUE)).collect(Collectors.toList()));
        }

        for (final String caller : SomaticVariantConstants.ALL_CALLERS) {
            checks.add(new HealthCheck(runContext.tumorSample(), SomaticCheck.AF_LOWER_SD.checkName(caller),
                    HealthCheckConstants.ERROR_VALUE));
            checks.add(new HealthCheck(runContext.tumorSample(), SomaticCheck.AF_MEDIAN.checkName(caller),
                    HealthCheckConstants.ERROR_VALUE));
            checks.add(new HealthCheck(runContext.tumorSample(), SomaticCheck.AF_UPPER_SD.checkName(caller),
                    HealthCheckConstants.ERROR_VALUE));
        }

        return toMultiValueResult(checks);
    }

    @NotNull
    private BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);
        return new MultiValueResult(checkType(), checks);
    }

    @NotNull
    private static List<SomaticVariant> toSomaticVariants(@NotNull final List<String> lines) {
        return lines.stream().map(SomaticVariantFactory::fromVCFLine).collect(Collectors.toList());
    }

    @NotNull
    private static List<HealthCheck> getTypeChecks(@NotNull final List<SomaticVariant> variants,
            @NotNull final String sampleId, @NotNull final VariantType type) {
        final List<HealthCheck> checks = new ArrayList<>();
        final List<SomaticVariant> variantsForType = VariantFilter.filter(variants, VariantPredicates.withType(type));

        final HealthCheck variantCountCheck = new HealthCheck(sampleId,
                SomaticCheck.COUNT_TOTAL.checkName(type.name()), String.valueOf(variantsForType.size()));
        checks.add(variantCountCheck);

        final List<SomaticVariant> variantsWithDBSNPAndNotCOSMIC = VariantFilter.filter(variantsForType,
                VariantPredicates.inDBSNPAndNotInCOSMIC());
        final HealthCheck dbsnpCheck = new HealthCheck(sampleId, SomaticCheck.DBSNP_COUNT.checkName(type.name()),
                String.valueOf(variantsWithDBSNPAndNotCOSMIC.size()));
        checks.add(dbsnpCheck);

        for (final String caller : SomaticVariantConstants.ALL_CALLERS) {
            List<SomaticVariant> callsPerCaller = VariantFilter.filter(variantsForType,
                    VariantPredicates.withCaller(caller));
            checks.add(new HealthCheck(sampleId, SomaticCheck.COUNT_PER_CALLER.checkName(type.name(), caller),
                    String.valueOf(callsPerCaller.size())));
        }

        final List<HealthCheck> precisionChecks = SomaticVariantConstants.ALL_CALLERS.stream().map(
                caller -> calculatePrecision(variantsForType, sampleId, type, caller)).collect(Collectors.toList());
        checks.addAll(precisionChecks);

        final List<HealthCheck> sensitivityChecks = SomaticVariantConstants.ALL_CALLERS.stream().map(
                caller -> calculateSensitivity(variantsForType, sampleId, type, caller)).collect(Collectors.toList());
        checks.addAll(sensitivityChecks);

        final List<HealthCheck> proportionChecks = CALLERS_COUNT.stream().map(
                callerCount -> calculateProportion(variantsForType, sampleId, type, callerCount)).collect(
                Collectors.toList());
        checks.addAll(proportionChecks);
        return checks;
    }

    @NotNull
    private static List<HealthCheck> getAFChecks(final List<SomaticVariant> variants, final String sampleId) {
        final List<HealthCheck> checks = Lists.newArrayList();
        for (String caller : SomaticVariantConstants.ALL_CALLERS) {
            List<SomaticVariant> filteredVariants = VariantFilter.filter(variants,
                    VariantPredicates.withCaller(caller));

            if (filteredVariants.size() > 0) {
                List<Double> alleleFreqs = filteredVariants.stream().map(SomaticVariant::alleleFrequency).collect(
                        Collectors.toList());
                alleleFreqs.sort(Comparator.naturalOrder());

                int lowerSDIndex = (int) Math.round(alleleFreqs.size() * AF_SD_DISTANCE);
                int medianIndex = Math.min(alleleFreqs.size() - 1, (int) Math.round(alleleFreqs.size() / 2D));
                int upperSDIndex = Math.min(alleleFreqs.size() - 1,
                        (int) Math.round(alleleFreqs.size() * (1 - AF_SD_DISTANCE)));

                checks.add(new HealthCheck(sampleId, SomaticCheck.AF_LOWER_SD.checkName(caller),
                        String.valueOf(alleleFreqs.get(lowerSDIndex))));
                checks.add(new HealthCheck(sampleId, SomaticCheck.AF_MEDIAN.checkName(caller),
                        String.valueOf(alleleFreqs.get(medianIndex))));
                checks.add(new HealthCheck(sampleId, SomaticCheck.AF_UPPER_SD.checkName(caller),
                        String.valueOf(alleleFreqs.get(upperSDIndex))));
            } else {
                checks.add(new HealthCheck(sampleId, SomaticCheck.AF_LOWER_SD.checkName(caller), "-"));
                checks.add(new HealthCheck(sampleId, SomaticCheck.AF_MEDIAN.checkName(caller), "-"));
                checks.add(new HealthCheck(sampleId, SomaticCheck.AF_UPPER_SD.checkName(caller), "-"));
            }
        }
        return checks;
    }

    @NotNull
    private static HealthCheck calculatePrecision(@NotNull final List<SomaticVariant> variants,
            @NotNull final String sampleId, @NotNull final VariantType type, @NotNull final String caller) {
        final List<SomaticVariant> variantsForCaller = VariantFilter.filter(variants,
                VariantPredicates.withCaller(caller));
        final List<SomaticVariant> variantsForCallerWithMoreThanOneCaller = VariantFilter.filter(variantsForCaller,
                VariantPredicates.withMinCallers(2));

        double precision = 0D;
        if (!variantsForCallerWithMoreThanOneCaller.isEmpty() && !variantsForCaller.isEmpty()) {
            precision = (double) variantsForCallerWithMoreThanOneCaller.size() / variantsForCaller.size();
        }
        return new HealthCheck(sampleId, SomaticCheck.PRECISION_CHECK.checkName(type.name(), caller),
                String.valueOf(precision));
    }

    @NotNull
    private static HealthCheck calculateSensitivity(@NotNull final List<SomaticVariant> variants,
            @NotNull final String sampleId, @NotNull final VariantType type, @NotNull final String caller) {
        final List<SomaticVariant> variantsWithMoreThanOneCaller = VariantFilter.filter(variants,
                VariantPredicates.withMinCallers(2));
        final List<SomaticVariant> variantsForCallerWithMoreThanOneCaller = VariantFilter.filter(
                variantsWithMoreThanOneCaller, VariantPredicates.withCaller(caller));

        double sensitivity = 0D;
        if (!variantsForCallerWithMoreThanOneCaller.isEmpty() && !variantsWithMoreThanOneCaller.isEmpty()) {
            sensitivity =
                    (double) variantsForCallerWithMoreThanOneCaller.size() / variantsWithMoreThanOneCaller.size();
        }
        return new HealthCheck(sampleId, SomaticCheck.SENSITIVITY_CHECK.checkName(type.name(), caller),
                String.valueOf(sensitivity));
    }

    @NotNull
    private static HealthCheck calculateProportion(@NotNull final List<SomaticVariant> variants,
            @NotNull final String sampleId, @NotNull final VariantType type, final int count) {
        final List<SomaticVariant> variantsWithCallerCount = VariantFilter.filter(variants,
                VariantPredicates.withExactCallerCount(count));
        double proportion = 0D;
        if (!variantsWithCallerCount.isEmpty() && !variants.isEmpty()) {
            proportion = (double) variantsWithCallerCount.size() / variants.size();
        }

        return new HealthCheck(sampleId, SomaticCheck.PROPORTION_CHECK.checkName(type.name(), String.valueOf(count)),
                String.valueOf(proportion));
    }
}
