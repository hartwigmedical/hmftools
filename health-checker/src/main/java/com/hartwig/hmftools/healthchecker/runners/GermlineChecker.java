package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.GermlineCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.GERMLINE)
public class GermlineChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(GermlineChecker.class);

    private static final String GERMLINE_VCF_EXTENSION_V1_9 = "_GoNLv5.vcf";
    private static final String GERMLINE_VCF_EXTENSION_V1_10 = ".annotated.vcf";

    public GermlineChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.GERMLINE;
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HartwigException {
        List<GermlineVariant> variants;
        try {
            variants = VCFFileLoader.loadGermlineVCF(runContext.runDirectory(), GERMLINE_VCF_EXTENSION_V1_10);
        } catch (IOException exception) {
            variants = VCFFileLoader.loadGermlineVCF(runContext.runDirectory(), GERMLINE_VCF_EXTENSION_V1_9);
        }

        variants = VariantFilter.passOnly(variants);
        final List<HealthCheck> refChecks = calcChecksForSample(variants, runContext.refSample(), true);
        if (runContext.isSomaticRun()) {
            final List<HealthCheck> tumorChecks = calcChecksForSample(variants, runContext.tumorSample(), false);

            return toPatientResult(refChecks, tumorChecks);
        } else {
            return toMultiValueResult(refChecks);
        }
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        if (runContext.isSomaticRun()) {
            return toPatientResult(getErrorChecksForSample(runContext.refSample()),
                    getErrorChecksForSample(runContext.tumorSample()));
        } else {
            return toMultiValueResult(getErrorChecksForSample(runContext.refSample()));
        }
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
    private BaseResult toPatientResult(@NotNull final List<HealthCheck> refChecks,
            @NotNull final List<HealthCheck> tumorChecks) {
        HealthCheck.log(LOGGER, refChecks);
        HealthCheck.log(LOGGER, tumorChecks);

        return new PatientResult(checkType(), refChecks, tumorChecks);
    }

    @NotNull
    private BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);

        return new MultiValueResult(checkType(), checks);
    }

    @NotNull
    private static List<HealthCheck> calcChecksForSample(@NotNull List<GermlineVariant> variants,
            @NotNull final String sampleId, final boolean isRefSample) {
        final HealthCheck snp = countValidVariantsPerType(sampleId, variants, VariantType.SNP,
                GermlineCheck.VARIANTS_GERMLINE_SNP, isRefSample);
        final HealthCheck indel = countValidVariantsPerType(sampleId, variants, VariantType.INDEL,
                GermlineCheck.VARIANTS_GERMLINE_INDELS, isRefSample);
        return Arrays.asList(snp, indel);
    }

    @NotNull
    private static HealthCheck countValidVariantsPerType(@NotNull final String sampleId,
            @NotNull final List<GermlineVariant> variants, @NotNull final VariantType type,
            @NotNull final GermlineCheck check, final boolean isRefSample) {
        final long count = VariantFilter.validGermline(variants, type, isRefSample).size();
        return new HealthCheck(sampleId, check.toString(), String.valueOf(count));
    }
}
