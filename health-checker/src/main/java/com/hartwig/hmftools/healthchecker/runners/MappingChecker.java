package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.flagstatreader.FlagStatData;
import com.hartwig.hmftools.healthchecker.flagstatreader.FlagStats;
import com.hartwig.hmftools.healthchecker.flagstatreader.FlagStatsType;
import com.hartwig.hmftools.healthchecker.flagstatreader.SambambaFlagStatParser;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;
import com.hartwig.hmftools.healthchecker.runners.checks.MappingCheck;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("WeakerAccess")
@ResourceWrapper(type = CheckType.MAPPING)
public class MappingChecker extends ErrorHandlingChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(MappingChecker.class);

    private static final String FASTQC_BASE_DIRECTORY = "QCStats";
    private static final String FLAGSTAT_BASE_DIRECTORY = "mapping";

    private static final String FLAGSTAT_FILE_FILTER = ".realign";

    public MappingChecker() {
    }

    @NotNull
    @Override
    public CheckType checkType() {
        return CheckType.MAPPING;
    }

    @NotNull
    @Override
    public BaseResult tryRun(@NotNull final RunContext runContext) throws IOException, HartwigException {
        final List<HealthCheck> refChecks = getChecksForSample(runContext.runDirectory(), runContext.refSample());
        final List<HealthCheck> tumorChecks = getChecksForSample(runContext.runDirectory(), runContext.tumorSample());

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
        final List<HealthCheck> checks = Lists.newArrayList();
        for (MappingCheck check : MappingCheck.values()) {
            checks.add(new HealthCheck(sampleId, check.toString(), HealthCheckConstants.ERROR_VALUE));
        }
        return checks;
    }

    @NotNull
    private BaseResult toPatientResult(@NotNull final List<HealthCheck> refChecks,
            @NotNull final List<HealthCheck> tumorChecks) {
        HealthCheck.log(LOGGER, refChecks);
        HealthCheck.log(LOGGER, tumorChecks);

        return new PatientResult(checkType(), refChecks, tumorChecks);
    }

    @NotNull
    private static List<HealthCheck> getChecksForSample(@NotNull final String runDirectory,
            @NotNull final String sampleId) throws IOException, HartwigException {
        final String basePathForTotalSequences = getBasePathForTotalSequences(runDirectory, sampleId);
        final long totalSequences = HealthCheckFunctions.sumOfTotalSequencesFromFastQC(basePathForTotalSequences);

        final String basePathForFlagStat = getBasePathForFlagStat(runDirectory, sampleId);
        return extractChecksFromFlagstats(basePathForFlagStat, sampleId, totalSequences);
    }

    @NotNull
    private static List<HealthCheck> extractChecksFromFlagstats(@NotNull final String basePath,
            @NotNull final String sampleId, final long totalSequences) throws IOException, EmptyFileException {
        final FlagStatData flagstatData = new SambambaFlagStatParser().parse(basePath, FLAGSTAT_FILE_FILTER);

        final List<FlagStats> passed = flagstatData.getPassedStats();

        final HealthCheck isAllReadDataReport = proportionOfReadsUsed(sampleId, totalSequences, passed);
        final HealthCheck mappedDataReport = mappedPercentage(sampleId, passed);
        final HealthCheck properlyPairedDataReport = properlyPairedPercentage(sampleId, passed);
        final HealthCheck singletonDataReport = singletonPercentage(sampleId, passed);
        final HealthCheck mateMappedDataReport = mateMappedDiffChrPercentage(sampleId, passed);
        final HealthCheck duplicateDataReport = duplicatePercentage(sampleId, passed);

        return Lists.newArrayList(isAllReadDataReport, mappedDataReport, properlyPairedDataReport, singletonDataReport,
                mateMappedDataReport, duplicateDataReport);
    }

    @NotNull
    private static HealthCheck mappedPercentage(@NotNull final String sampleId,
            @NotNull final List<FlagStats> passed) {
        final FlagStats mappedStat = passed.get(FlagStatsType.MAPPED_INDEX.getIndex());
        final FlagStats totalStat = passed.get(FlagStatsType.TOTAL_INDEX.getIndex());
        final double mappedPercentage = mappedStat.getValue() / totalStat.getValue();

        return new HealthCheck(sampleId, MappingCheck.MAPPING_PERCENTAGE_MAPPED.toString(),
                String.valueOf(mappedPercentage));
    }

    @NotNull
    private static HealthCheck properlyPairedPercentage(@NotNull final String sampleId,
            @NotNull final List<FlagStats> passed) {
        final FlagStats mappedStat = passed.get(FlagStatsType.MAPPED_INDEX.getIndex());
        final FlagStats properPaired = passed.get(FlagStatsType.PROPERLY_PAIRED_INDEX.getIndex());
        final double properlyPairedPercentage = properPaired.getValue() / mappedStat.getValue();

        return new HealthCheck(sampleId, MappingCheck.MAPPING_PROPERLY_PAIRED_PROPORTION_OF_MAPPED.toString(),
                String.valueOf(properlyPairedPercentage));
    }

    @NotNull
    private static HealthCheck singletonPercentage(@NotNull final String sampleId,
            @NotNull final List<FlagStats> passed) {
        final FlagStats mappedStat = passed.get(FlagStatsType.MAPPED_INDEX.getIndex());
        final FlagStats singletonStat = passed.get(FlagStatsType.SINGLETONS_INDEX.getIndex());
        final double singletonPercentage = singletonStat.getValue() / mappedStat.getValue();

        return new HealthCheck(sampleId, MappingCheck.MAPPING_PROPORTION_SINGLETON.toString(),
                String.valueOf(singletonPercentage));
    }

    @NotNull
    private static HealthCheck mateMappedDiffChrPercentage(@NotNull final String sampleId,
            @NotNull final List<FlagStats> passed) {
        final FlagStats mappedStat = passed.get(FlagStatsType.MAPPED_INDEX.getIndex());
        final FlagStats diffPercStat = passed.get(FlagStatsType.MATE_MAP_DIF_CHR_INDEX.getIndex());
        final double mateMappedDiffChrPerc = diffPercStat.getValue() / mappedStat.getValue();

        return new HealthCheck(sampleId, MappingCheck.MAPPING_PROPORTION_MAPPED_DIFFERENT_CHR.toString(),
                String.valueOf(mateMappedDiffChrPerc));
    }

    @NotNull
    private static HealthCheck duplicatePercentage(@NotNull final String sampleId,
            @NotNull final List<FlagStats> passed) {
        final FlagStats totalStat = passed.get(FlagStatsType.TOTAL_INDEX.getIndex());
        final FlagStats duplicateStat = passed.get(FlagStatsType.DUPLICATES_INDEX.getIndex());
        final double proportionOfDuplicateRead = duplicateStat.getValue() / totalStat.getValue();

        return new HealthCheck(sampleId, MappingCheck.MAPPING_MARKDUP_PROPORTION_DUPLICATES.toString(),
                String.valueOf(proportionOfDuplicateRead));
    }

    @NotNull
    private static HealthCheck proportionOfReadsUsed(@NotNull final String sampleId, final long totalSequences,
            @NotNull final List<FlagStats> passed) {
        final FlagStats totalStat = passed.get(FlagStatsType.TOTAL_INDEX.getIndex());
        final FlagStats secondaryStat = passed.get(FlagStatsType.SECONDARY_INDEX.getIndex());
        double proportionReadPercentage = (totalStat.getValue() - secondaryStat.getValue()) / totalSequences;

        return new HealthCheck(sampleId, MappingCheck.MAPPING_PROPORTION_READ_VS_TOTAL_SEQUENCES.toString(),
                String.valueOf(proportionReadPercentage));
    }

    @NotNull
    private static String getBasePathForFlagStat(@NotNull final String runDirectory, @NotNull final String sampleId) {
        return runDirectory + File.separator + sampleId + File.separator + FLAGSTAT_BASE_DIRECTORY;
    }

    @NotNull
    private static String getBasePathForTotalSequences(@NotNull final String runDirectory,
            @NotNull final String sampleId) {
        return runDirectory + File.separator + sampleId + File.separator + FASTQC_BASE_DIRECTORY;
    }
}
