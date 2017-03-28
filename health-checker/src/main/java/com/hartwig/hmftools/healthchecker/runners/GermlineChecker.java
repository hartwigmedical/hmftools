package com.hartwig.hmftools.healthchecker.runners;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcf.VCFFileStreamer;
import com.hartwig.hmftools.healthchecker.context.RunContext;
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
        final BufferedReader reader = openReader(runContext.runDirectory());

        final Map<VariantType, Long> refSampleCountPerType = buildInitialCountMap();
        final Map<VariantType, Long> tumorSampleCountPerType = buildInitialCountMap();

        GermlineVariant variant = VCFFileStreamer.nextVariant(reader);
        while (variant != null) {
            if (VariantFilter.isPass(variant)) {
                for (final VariantType type : typesToInclude()) {
                    if (VariantFilter.validGermline(variant, type, true)) {
                        refSampleCountPerType.put(type, refSampleCountPerType.get(type) + 1L);
                    }
                    if (runContext.isSomaticRun() && VariantFilter.validGermline(variant, type, false)) {
                        tumorSampleCountPerType.put(type, tumorSampleCountPerType.get(type) + 1L);
                    }
                }
            }
            variant = VCFFileStreamer.nextVariant(reader);
        }

        final List<HealthCheck> refChecks = buildChecks(runContext.refSample(), refSampleCountPerType);

        if (runContext.isSomaticRun()) {
            final List<HealthCheck> tumorChecks = buildChecks(runContext.tumorSample(), tumorSampleCountPerType);
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
    private static BufferedReader openReader(@NotNull final String runDirectory) throws IOException, HartwigException {
        try {
            return VCFFileStreamer.getVCFReader(runDirectory, GERMLINE_VCF_EXTENSION_V1_10);
        } catch (FileNotFoundException exception) {
            return VCFFileStreamer.getVCFReader(runDirectory, GERMLINE_VCF_EXTENSION_V1_9);
        }
    }

    @NotNull
    private static Map<VariantType, Long> buildInitialCountMap() {
        final Map<VariantType, Long> initialCountMap = Maps.newHashMap();
        for (final VariantType type : typesToInclude()) {
            initialCountMap.put(type, 0L);
        }
        return initialCountMap;
    }

    @NotNull
    private static List<HealthCheck> buildChecks(@NotNull final String sample,
            @NotNull final Map<VariantType, Long> countsPerType) {
        return Lists.newArrayList(new HealthCheck(sample, GermlineCheck.VARIANTS_GERMLINE_SNP.toString(),
                        Long.toString(countsPerType.get(VariantType.SNP))),
                new HealthCheck(sample, GermlineCheck.VARIANTS_GERMLINE_INDELS.toString(),
                        Long.toString(countsPerType.get(VariantType.INDEL))));
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
    private static List<VariantType> typesToInclude() {
        return Lists.newArrayList(VariantType.SNP, VariantType.INDEL);
    }
}
