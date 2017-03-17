package com.hartwig.hmftools.patientreporter.batch;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.patientreporter.algo.GenomeAnalysis;
import com.hartwig.hmftools.patientreporter.algo.SinglePatientReporter;
import com.hartwig.hmftools.patientreporter.util.FindingsToCSV;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BatchReportAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(SinglePatientReporter.class);

    @NotNull
    private final SinglePatientReporter singlePatientReporter;
    @Nullable
    private final String statsFile;
    @Nullable
    private final String variantFindingsFile;
    @Nullable
    private final String copyNumberFindingsFile;

    public BatchReportAnalyser(@NotNull final SinglePatientReporter singlePatientReporter,
            @Nullable final String statsFile, @Nullable final String variantFindingsFile,
            @Nullable final String copyNumberFindingsFile) {
        this.singlePatientReporter = singlePatientReporter;
        this.statsFile = statsFile;
        this.variantFindingsFile = variantFindingsFile;
        this.copyNumberFindingsFile = copyNumberFindingsFile;
    }

    @NotNull
    public Map<String, RunStats> run(@NotNull final String batchDirectory) throws IOException, HartwigException {
        final Map<String, RunStats> runStatsPerSample = generateStatsForBatch(batchDirectory);

        if (statsFile != null) {
            LOGGER.info("Writing stats to " + statsFile);
            writeStatsToOutputFile(runStatsPerSample, statsFile);
        }

        if (variantFindingsFile != null) {
            LOGGER.info("Writing findings to " + variantFindingsFile);
            writeVariantFindings(runStatsPerSample, variantFindingsFile);
        }

        if (copyNumberFindingsFile != null) {
            LOGGER.info("Writing findings to " + copyNumberFindingsFile);
            writeCopyNumberFindings(runStatsPerSample, copyNumberFindingsFile);
        }

        return runStatsPerSample;
    }

    @NotNull
    private Map<String, RunStats> generateStatsForBatch(@NotNull final String directory)
            throws IOException, HartwigException {
        final Map<String, RunStats> statsPerSample = Maps.newHashMap();

        LOGGER.info("Running patient reporter batch-mode on " + directory);
        for (final Path run : Files.list(new File(directory).toPath()).collect(Collectors.toList())) {
            LOGGER.info(" Running on " + run.toFile().getPath());
            final GenomeAnalysis genomeAnalysis = singlePatientReporter.analyseGenomeData(run.toFile().getPath());

            statsPerSample.put(genomeAnalysis.sample(), RunStats.fromGenomeAnalysis(genomeAnalysis));
        }
        return statsPerSample;
    }

    private static void writeStatsToOutputFile(@NotNull final Map<String, RunStats> statsPerSample,
            @NotNull final String outputFile) throws IOException {
        final VariantConsequence[] consequences = VariantConsequence.values();
        final List<String> output = Lists.newArrayList();

        String header = "SAMPLE,VARIANT_COUNT,PASS_ONLY_COUNT,CONSENSUS_COUNT,MUTATIONAL_LOAD,CONSEQUENCE_COUNT,CNV_COUNT";
        for (final VariantConsequence consequence : consequences) {
            header += ("," + consequence.name() + "_COUNT");
        }
        output.add(header);

        for (final Map.Entry<String, RunStats> entry : statsPerSample.entrySet()) {
            final RunStats stats = entry.getValue();

            String consequenceList = Strings.EMPTY;
            for (final VariantConsequence consequence : consequences) {
                consequenceList += ("," + stats.consequenceCounts().get(consequence).toString());
            }

            output.add(entry.getKey() + "," + stats.allVariantCount() + "," + stats.passedVariantCount() + ","
                    + stats.consensusPassedCount() + "," + stats.mutationalLoad() + ","
                    + stats.variantFindings().size() + "," + stats.copyNumberFindings().size() + consequenceList);
        }
        Files.write(new File(outputFile).toPath(), output);
    }

    private static void writeVariantFindings(@NotNull final Map<String, RunStats> statsPerSample,
            @NotNull final String outputFile) throws IOException {
        final List<String> output = Lists.newArrayList();

        for (final Map.Entry<String, RunStats> entry : statsPerSample.entrySet()) {
            appendFindings(entry.getKey(), output, FindingsToCSV.varToCSV(entry.getValue().variantFindings()));
        }
        Files.write(new File(outputFile).toPath(), output);
    }

    private static void writeCopyNumberFindings(@NotNull final Map<String, RunStats> statsPerSample,
            @NotNull final String outputFile) throws IOException {
        final List<String> output = Lists.newArrayList();

        for (final Map.Entry<String, RunStats> entry : statsPerSample.entrySet()) {
            appendFindings(entry.getKey(), output, FindingsToCSV.cnvToCSV(entry.getValue().copyNumberFindings()));
        }

        Files.write(new File(outputFile).toPath(), output);
    }

    @NotNull
    private static List<String> appendFindings(@NotNull final String sample, @NotNull final List<String> current,
            @NotNull final List<String> findings) {
        if (current.isEmpty()) {
            current.add("SAMPLE," + findings.get(0));

        }
        current.addAll(findings.subList(1, findings.size()).stream().map(finding -> sample + "," + finding).collect(
                Collectors.toList()));
        return current;
    }
}
