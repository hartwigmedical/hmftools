package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.copynumber.CopyNumberFactory;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathPrefixSuffixFinder;
import com.hartwig.hmftools.common.io.reader.FileReader;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.runners.checks.CopynumberCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.COPYNUMBER)
public class CopynumberChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(CopynumberChecker.class);

    // KODU: copynumber data is stored in {run}/copyNumber/{sampleR}_{sampleT}/freec/{sampleT}<>.bam_CNVs
    private static final String COPYNUMBER_BASE_DIRECTORY = "copyNumber";
    private static final String COPYNUMBER_SAMPLE_CONNECTOR = "_";
    private static final String COPYNUMBER_ALGO_DIRECTORY = "freec";
    private static final String COPYNUMBER_SUFFIX = ".bam_CNVs";
    private static final String COPYNUMBER_RATIO_SUFFIX = ".bam_ratio.txt";

    public CopynumberChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.COPYNUMBER;
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HartwigException {
        long totalGain = 0;
        long totalLoss = 0;
        for (final CopyNumber copyNumber : copynumberLines(runContext)) {
            if (copyNumber.isGain()) {
                totalGain += copyNumber.basesAffected();
            } else if (copyNumber.isLoss()) {
                totalLoss += copyNumber.basesAffected();
            }
        }

        return toMultiValueResult(runContext, String.valueOf(totalGain), String.valueOf(totalLoss));
    }

    @NotNull
    private static List<CopyNumber> copynumberLines(@NotNull final RunContext runContext)
            throws IOException, HartwigException {
        final Path copynumberPath = PathPrefixSuffixFinder.build().findPath(getBasePath(runContext),
                relevantSample(runContext), COPYNUMBER_SUFFIX);
        final List<String> lines;
        try {
            lines = FileReader.build().readLines(copynumberPath);
        } catch (EmptyFileException e) {
            // if the CNV is empty (but exists) and the ratio file exists, there is no problem (just no CNVs found)
            final Path ratioPath = PathPrefixSuffixFinder.build().findPath(getBasePath(runContext),
                    relevantSample(runContext), COPYNUMBER_RATIO_SUFFIX);
            FileReader.build().readLines(ratioPath);
            return Collections.emptyList();
        }

        List<CopyNumber> copyNumbers = Lists.newArrayList();
        for (String line : lines) {
            copyNumbers.add(CopyNumberFactory.fromCNVLine(line));
        }
        return copyNumbers;
    }

    @NotNull
    @Override
    public BaseResult errorRun(@NotNull final RunContext runContext) {
        return toMultiValueResult(runContext, HealthCheckConstants.ERROR_VALUE, HealthCheckConstants.ERROR_VALUE);
    }

    @NotNull
    private BaseResult toMultiValueResult(@NotNull final RunContext runContext, @NotNull final String totalGain,
            @NotNull final String totalLoss) {
        final HealthCheck gainCheck = new HealthCheck(relevantSample(runContext),
                CopynumberCheck.COPYNUMBER_GENOME_GAIN.toString(), totalGain);
        final HealthCheck lossCheck = new HealthCheck(relevantSample(runContext),
                CopynumberCheck.COPYNUMBER_GENOME_LOSS.toString(), totalLoss);
        final List<HealthCheck> checks = Lists.newArrayList(gainCheck, lossCheck);
        HealthCheck.log(LOGGER, checks);
        return new MultiValueResult(checkType(), checks);
    }

    @NotNull
    private static String getBasePath(@NotNull final RunContext runContext) {
        final String baseDir = runContext.runDirectory() + File.separator + COPYNUMBER_BASE_DIRECTORY + File.separator;
        final String sampleDir = runContext.isSomaticRun() ?
                runContext.refSample() + COPYNUMBER_SAMPLE_CONNECTOR + runContext.tumorSample() :
                runContext.refSample();

        return baseDir + sampleDir + File.separator + COPYNUMBER_ALGO_DIRECTORY;
    }

    @NotNull
    private static String relevantSample(@NotNull final RunContext runContext) {
        return runContext.isSomaticRun() ? runContext.tumorSample() : runContext.refSample();
    }
}
