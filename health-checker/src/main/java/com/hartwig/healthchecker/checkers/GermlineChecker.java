package com.hartwig.healthchecker.checkers;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.healthchecker.checkers.checks.GermlineCheck;
import com.hartwig.healthchecker.checkers.checks.HealthCheck;
import com.hartwig.hmftools.common.variant.VCFGermlineData;
import com.hartwig.hmftools.common.variant.VCFGermlineDataFactory;
import com.hartwig.hmftools.common.variant.VCFType;
import com.hartwig.hmftools.common.variant.predicate.VCFGermlineVariantPredicate;
import com.hartwig.healthchecker.exception.HealthChecksException;
import com.hartwig.healthchecker.io.dir.RunContext;
import com.hartwig.healthchecker.io.path.PathExtensionFinder;
import com.hartwig.healthchecker.io.reader.LineReader;
import com.hartwig.hmftools.common.variant.predicate.VCFPassDataLinePredicate;
import com.hartwig.healthchecker.resource.ResourceWrapper;
import com.hartwig.healthchecker.result.BaseResult;
import com.hartwig.healthchecker.result.PatientResult;

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
        final List<VCFGermlineData> variants = getVCFDataForGermLine(passFilterLines);

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
    private static List<HealthCheck> calcChecksForSample(@NotNull List<VCFGermlineData> vcfData,
            @NotNull final String sampleId, final boolean isRefSample) {
        final HealthCheck snp = getGermlineVariantCount(sampleId, vcfData, VCFType.SNP,
                GermlineCheck.VARIANTS_GERMLINE_SNP, isRefSample);
        final HealthCheck indels = getGermlineVariantCount(sampleId, vcfData, VCFType.INDELS,
                GermlineCheck.VARIANTS_GERMLINE_INDELS, isRefSample);
        return Arrays.asList(snp, indels);
    }

    @NotNull
    private static List<VCFGermlineData> getVCFDataForGermLine(@NotNull final List<String> lines) {
        return lines.stream().map(VCFGermlineDataFactory::fromVCFLine).filter(notNull()).collect(Collectors.toList());
    }

    @NotNull
    private static HealthCheck getGermlineVariantCount(@NotNull final String sampleId,
            @NotNull final List<VCFGermlineData> vcfData, @NotNull final VCFType vcfType,
            @NotNull final GermlineCheck check, final boolean isRefSample) {
        final long count = vcfData.stream().filter(new VCFGermlineVariantPredicate(vcfType, isRefSample)).count();
        return new HealthCheck(sampleId, check.toString(), String.valueOf(count));
    }

    @NotNull
    private static Predicate<VCFGermlineData> notNull() {
        return vcf -> vcf != null;
    }
}
