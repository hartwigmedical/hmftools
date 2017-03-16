package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.patientreporter.util.ConsequenceCount;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class BatchReportAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(SinglePatientReporter.class);

    @NotNull
    private final SinglePatientReporter singlePatientReporter;

    BatchReportAnalyser(@NotNull final SinglePatientReporter singlePatientReporter) {
        this.singlePatientReporter = singlePatientReporter;
    }

    @NotNull
    List<String> run(@NotNull final String batchDirectory) throws IOException, HartwigException {
        final VariantConsequence[] consequences = VariantConsequence.values();
        final List<String> output = Lists.newArrayList();

        String header = "SAMPLE,VARIANT_COUNT,PASS_ONLY_COUNT,CONSENSUS_COUNT,MISSENSE_COUNT,CONSEQUENCE_COUNT";
        for (final VariantConsequence consequence : consequences) {
            header += ("," + consequence.name() + "_COUNT");
        }
        output.add(header);
        LOGGER.info("Running patient reporter batch-mode on " + batchDirectory);
        for (final Path run : Files.list(new File(batchDirectory).toPath()).collect(Collectors.toList())) {
            LOGGER.info(" Running on " + run.toFile().getPath());
            final VariantAnalysis analysis = singlePatientReporter.runVariantAnalysis(run.toFile().getPath());

            final Map<VariantConsequence, Integer> counts = ConsequenceCount.count(analysis.consensusPassedVariants());
            String consequenceList = Strings.EMPTY;
            for (final VariantConsequence consequence : consequences) {
                consequenceList += ("," + counts.get(consequence).toString());
            }

            output.add(analysis.sample() + "," + analysis.allVariants().size() + "," + analysis.passedVariants().size()
                    + "," + analysis.consensusPassedVariants().size() + "," + analysis.missenseVariants().size() + ","
                    + analysis.findings().size() + consequenceList);
        }

        return output;
    }
}
