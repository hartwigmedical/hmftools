package com.hartwig.hmftools.summon;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.summon.actionability.ActionabilityEntry;
import com.hartwig.hmftools.summon.actionability.ActionabilityFileReader;
import com.hartwig.hmftools.summon.conclusion.ActionabilityConclusion;
import com.hartwig.hmftools.summon.conclusion.ConclusionAlgo;
import com.hartwig.hmftools.summon.conclusion.SummonConclusionFile;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SummonApplication {

    private static final Logger LOGGER = LogManager.getLogger(SummonApplication.class);
    private static final String VERSION = SummonApplication.class.getPackage().getImplementationVersion();

    public static void main(@NotNull String[] args) throws IOException {
        LOGGER.info("Running SUMMON v{}", VERSION);

        Options options = SummonConfig.createOptions();

        SummonConfig config = null;
        try {
            config = SummonConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("SUMMON", options);
            System.exit(1);
        }

        new SummonApplication(config).run();

        LOGGER.info("Complete");
    }

    @NotNull
    private final SummonConfig config;

    private SummonApplication(@NotNull final SummonConfig config) {
        this.config = config;
    }

    public void run() throws IOException {
        LOGGER.info("Running SUMMON algo on sample {})", config.tumorSampleId());

        SummonAlgo algo = SummonAlgo.build(config.actionabilityDatabaseTsv(), config.driverGene37Tsv(), config.driverGene38Tsv());
        SummonData summonData = algo.run(config);

        ActionabilityConclusion actionabilityConclusion = ConclusionAlgo.generateConclusion(summonData);

        String filename = SummonConclusionFile.generateFilename(config.outputDir(), config.tumorSampleId());
        LOGGER.info("Writing actionability conclusion to file: {}", filename);
        SummonConclusionFile.write(filename, actionabilityConclusion);
    }
}