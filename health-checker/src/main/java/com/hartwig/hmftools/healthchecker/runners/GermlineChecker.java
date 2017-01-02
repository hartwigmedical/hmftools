package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HealthChecksException;
import com.hartwig.hmftools.common.io.dir.RunContext;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.io.reader.LineReader;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.predicate.VCFGermlineVariantPredicate;
import com.hartwig.hmftools.common.variant.predicate.VCFPassDataLinePredicate;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.GermlineCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.GERMLINE)
public class GermlineChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(GermlineChecker.class);

    private static final String GERMLINE_VCF_EXTENSION_V1_9 = "_GoNLv5.vcf";
    private static final String GERMLINE_VCF_EXTENSION_V1_10 = ".filtered_variants.annotated.vcf";

    public GermlineChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.GERMLINE;
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HealthChecksException {
        Path vcfPath;
        try {
            vcfPath = PathExtensionFinder.build().findPath(runContext.runDirectory(), GERMLINE_VCF_EXTENSION_V1_10);
        } catch (IOException exception) {
            vcfPath = PathExtensionFinder.build().findPath(runContext.runDirectory(), GERMLINE_VCF_EXTENSION_V1_9);
        }

        final List<String> passFilterLines = LineReader.build().readLines(vcfPath, new VCFPassDataLinePredicate());
        final List<GermlineVariant> variants = getVariantsForGermLine(passFilterLines);

        final List<HealthCheck> refChecks = calcChecksForSample(variants, runContext.refSample(), true);
        final List<HealthCheck> tumorChecks = calcChecksForSample(variants, runContext.tumorSample(), false);

        return toPatientResult(refChecks, tumorChecks);
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        return toPatientResult(getErrorChecksForSample(runContext.refSample()),
                getErrorChecksForSample(runContext.tumorSample()));
    }

    @NotNull
    private static List<HealthCheck> getErrorChecksForSample(@NotNull final String sampleId) {
        List<HealthCheck> errorChecks = Lists.newArrayList();
        for (GermlineCheck check : GermlineCheck.values()) {
            errorChecks.add(new HealthCheck(sampleId, check.toString(), HealthCheckConstants.ERROR_VALUE));
        }

        return errorChecks;
    }

    @NotNull
    private PatientResult toPatientResult(@NotNull final List<HealthCheck> refChecks,
            @NotNull final List<HealthCheck> tumorChecks) {
        HealthCheck.log(LOGGER, refChecks);
        HealthCheck.log(LOGGER, tumorChecks);

        return new PatientResult(checkType(), refChecks, tumorChecks);
    }

    @NotNull
    private static List<HealthCheck> calcChecksForSample(@NotNull List<GermlineVariant> vcfData,
            @NotNull final String sampleId, final boolean isRefSample) {
        final HealthCheck snp = getGermlineVariantCount(sampleId, vcfData, VariantType.SNP,
                GermlineCheck.VARIANTS_GERMLINE_SNP, isRefSample);
        final HealthCheck indel = getGermlineVariantCount(sampleId, vcfData, VariantType.INDEL,
                GermlineCheck.VARIANTS_GERMLINE_INDELS, isRefSample);
        return Arrays.asList(snp, indel);
    }

    @NotNull
    private static List<GermlineVariant> getVariantsForGermLine(@NotNull final List<String> lines) {
        return lines.stream().map(GermlineVariantFactory::fromVCFLine).filter(notNull()).collect(Collectors.toList());
    }

    @NotNull
    private static HealthCheck getGermlineVariantCount(@NotNull final String sampleId,
            @NotNull final List<GermlineVariant> variants, @NotNull final VariantType variantType,
            @NotNull final GermlineCheck check, final boolean isRefSample) {
        final long count = variants.stream().filter(new VCFGermlineVariantPredicate(variantType, isRefSample)).count();
        return new HealthCheck(sampleId, check.toString(), String.valueOf(count));
    }

    @NotNull
    private static Predicate<GermlineVariant> notNull() {
        return vcf -> vcf != null;
    }
}
