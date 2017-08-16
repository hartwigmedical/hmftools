package com.hartwig.hmftools.healthchecker.runners;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.GermlineSampleData;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcf.VCFFileStreamer;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
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
    private static final String GERMLINE_VCF_EXTENSION = ".annotated.vcf";

    private static final List<VariantType> TYPES_TO_INCLUDE = Lists.newArrayList(VariantType.SNP, VariantType.INDEL);
    private static final String HETEROZYGOUS_GENOTYPE = "0/1";

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
        final BufferedReader reader = openReader(runContext.runDirectory());

        final GermlineStats refStats = new GermlineStats();
        final GermlineStats tumorStats = new GermlineStats();

        GermlineVariant variant = VCFFileStreamer.nextVariant(reader);
        while (variant != null) {
            if (VariantFilter.isPass(variant) && TYPES_TO_INCLUDE.contains(variant.type())) {
                if (VariantFilter.validGermlineData(variant.refData())) {
                    updateStats(refStats, variant.type(), variant.refData());
                }

                GermlineSampleData tumorData = variant.tumorData();
                if (tumorData != null && VariantFilter.validGermlineData(tumorData)) {
                    updateStats(tumorStats, variant.type(), tumorData);
                }
            }

            variant = VCFFileStreamer.nextVariant(reader);
        }

        final List<HealthCheck> refChecks = buildChecks(runContext.refSample(), refStats);

        if (runContext.isSomaticRun()) {
            final List<HealthCheck> tumorChecks = buildChecks(runContext.tumorSample(), tumorStats);
            return toPatientResult(refChecks, tumorChecks);
        } else {
            return toMultiValueResult(refChecks);
        }
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        if (runContext.isSomaticRun()) {
            return toPatientResult(getErrorChecksForSample(runContext.refSample()), getErrorChecksForSample(runContext.tumorSample()));
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
    private static BufferedReader openReader(@NotNull final String runDirectory) throws IOException, HartwigException {
        return VCFFileStreamer.getVCFReader(runDirectory, GERMLINE_VCF_EXTENSION);
    }

    private static void updateStats(@NotNull final GermlineStats stats, @NotNull final VariantType type,
            @NotNull final GermlineSampleData sampleData) {
        assert TYPES_TO_INCLUDE.contains(type);

        if (type.equals(VariantType.SNP)) {
            stats.snpCount++;
        } else if (type.equals(VariantType.INDEL)) {
            stats.indelCount++;
        }

        if (sampleData.genoType().equals(HETEROZYGOUS_GENOTYPE)) {
            stats.heterozygousCount++;
            if (sampleData.alleleFrequency() > 0.5) {
                stats.heterozygousCountAbove50VAF++;
            } else if (sampleData.alleleFrequency() < 0.5) {
                stats.heterozygousCountBelow50VAF++;
            }
        }
    }

    @NotNull
    private static List<HealthCheck> buildChecks(@NotNull final String sample, @NotNull final GermlineStats stats) {
        return Lists.newArrayList(new HealthCheck(sample, GermlineCheck.GERMLINE_SNP_COUNT.toString(), Long.toString(stats.snpCount)),
                new HealthCheck(sample, GermlineCheck.GERMLINE_INDEL_COUNT.toString(), Long.toString(stats.indelCount)),
                new HealthCheck(sample, GermlineCheck.GERMLINE_HETEROZYGOUS_COUNT.toString(), Long.toString(stats.heterozygousCount)),
                new HealthCheck(sample, GermlineCheck.GERMLINE_HETEROZYGOUS_COUNT_ABOVE_50_VAF.toString(),
                        Long.toString(stats.heterozygousCountAbove50VAF)),
                new HealthCheck(sample, GermlineCheck.GERMLINE_HETEROZYGOUS_COUNT_BELOW_50_VAF.toString(),
                        Long.toString(stats.heterozygousCountBelow50VAF)));
    }

    @NotNull
    private BaseResult toPatientResult(@NotNull final List<HealthCheck> refChecks, @NotNull final List<HealthCheck> tumorChecks) {
        HealthCheck.log(LOGGER, refChecks);
        HealthCheck.log(LOGGER, tumorChecks);

        return new PatientResult(checkType(), refChecks, tumorChecks);
    }

    @NotNull
    private BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);

        return new MultiValueResult(checkType(), checks);
    }

    private static class GermlineStats {
        private long snpCount;
        private long indelCount;
        private long heterozygousCount;
        private long heterozygousCountAbove50VAF;
        private long heterozygousCountBelow50VAF;
    }
}