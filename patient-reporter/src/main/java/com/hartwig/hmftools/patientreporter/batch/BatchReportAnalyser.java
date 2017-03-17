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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class BatchReportAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(SinglePatientReporter.class);

    @NotNull
    private final SinglePatientReporter singlePatientReporter;

    public BatchReportAnalyser(@NotNull final SinglePatientReporter singlePatientReporter) {
        this.singlePatientReporter = singlePatientReporter;
    }

    @NotNull
    public List<String> run(@NotNull final String batchDirectory) throws IOException, HartwigException {
        final VariantConsequence[] consequences = VariantConsequence.values();
        final List<String> output = Lists.newArrayList();

        String header = "SAMPLE,VARIANT_COUNT,PASS_ONLY_COUNT,CONSENSUS_COUNT,MUTATIONAL_LOAD,CONSEQUENCE_COUNT";
        for (final VariantConsequence consequence : consequences) {
            header += ("," + consequence.name() + "_COUNT");
        }
        output.add(header);

        final Map<String, RunStats> runStatsPerSample = generateStatsForBatch(batchDirectory);
        for (final Map.Entry<String, RunStats> entry : runStatsPerSample.entrySet()) {
            final RunStats stats = entry.getValue();

            String consequenceList = Strings.EMPTY;
            for (final VariantConsequence consequence : consequences) {
                consequenceList += ("," + stats.consequenceCounts().get(consequence).toString());
            }

            output.add(entry.getKey() + "," + stats.allVariantCount() + "," + stats.passedVariantCount() + ","
                    + stats.consensusPassedCount() + "," + stats.mutationalLoad() + ","
                    + stats.variantFindings().size() + consequenceList);
        }

        return output;
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
}
