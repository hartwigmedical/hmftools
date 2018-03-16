package com.hartwig.hmftools.healthchecker;

import java.io.IOException;
import java.util.Collection;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.io.FolderChecker;
import com.hartwig.hmftools.healthchecker.report.JsonReport;
import com.hartwig.hmftools.healthchecker.report.Report;
import com.hartwig.hmftools.healthchecker.runners.AmberChecker;
import com.hartwig.hmftools.healthchecker.runners.CoverageChecker;
import com.hartwig.hmftools.healthchecker.runners.HealthChecker;
import com.hartwig.hmftools.healthchecker.runners.KinshipChecker;
import com.hartwig.hmftools.healthchecker.runners.PurpleChecker;
import com.hartwig.hmftools.healthchecker.runners.SomaticVariantsChecker;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HealthChecksApplication {

    private static final Logger LOGGER = LogManager.getLogger(HealthChecksApplication.class);

    private static final String RUN_DIRECTORY = "run_dir";
    private static final String REPORT_FILE_PATH = "report_file_path";

    @NotNull
    private final RunContext runContext;
    @NotNull
    private final String reportFilePath;

    private HealthChecksApplication(@NotNull final RunContext runContext, @NotNull final String reportFilePath) {
        this.runContext = runContext;
        this.reportFilePath = reportFilePath;
    }

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        final String reportFilePath = cmd.getOptionValue(REPORT_FILE_PATH);

        if (runDirectory == null || reportFilePath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks", options);
            System.exit(1);
        }

        runDirectory = FolderChecker.build().checkFolder(runDirectory);
        RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDirectory);

        new HealthChecksApplication(runContext, reportFilePath).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUN_DIRECTORY, true, "The path containing the data for a single run");
        options.addOption(REPORT_FILE_PATH, true, "The path where the report will be written to.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private void run() throws IOException {
        final Report report = new JsonReport();
        final Collection<HealthChecker> checkers = Lists.newArrayList(new CoverageChecker(),
                new SomaticVariantsChecker(),
                new KinshipChecker(),
                new PurpleChecker(),
                new AmberChecker());

        for (final HealthChecker checker : checkers) {
            report.addResult(checker.run(runContext));
        }

        final Optional<String> reportPath = report.generateReport(reportFilePath);
        reportPath.ifPresent(path -> LOGGER.info(String.format("Report generated -> \n%s", path)));
    }
}
