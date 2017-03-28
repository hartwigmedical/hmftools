package com.hartwig.hmftools.healthchecker;

import java.io.IOException;
import java.util.Collection;
import java.util.Optional;

import com.hartwig.hmftools.common.exception.GenerateReportException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.healthchecker.context.ProductionRunContextFactory;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.io.FolderChecker;
import com.hartwig.hmftools.healthchecker.report.HealthCheckReportFactory;
import com.hartwig.hmftools.healthchecker.report.Report;
import com.hartwig.hmftools.healthchecker.runners.HealthChecker;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import rx.Observable;
import rx.functions.Action1;
import rx.observables.BlockingObservable;
import rx.schedulers.Schedulers;

public final class HealthChecksApplication {

    private static final Logger LOGGER = LogManager.getLogger(HealthChecksApplication.class);

    private static final String RUN_DIRECTORY = "run_dir";
    private static final String REPORT_TYPE = "report_type";
    private static final String REPORT_OUTPUT_PATH = "report_output_path";

    @NotNull
    private final RunContext runContext;
    @NotNull
    private final String reportType;
    @NotNull
    private final String reportOutputPath;

    private HealthChecksApplication(@NotNull final RunContext runContext, @NotNull final String reportType,
            @NotNull final String reportOutputPath) {
        this.runContext = runContext;
        this.reportType = reportType;
        this.reportOutputPath = reportOutputPath;
    }

    public static void main(final String... args) throws ParseException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        final String reportType = cmd.getOptionValue(REPORT_TYPE);
        final String reportOutputPath = cmd.getOptionValue(REPORT_OUTPUT_PATH);

        if (runDirectory == null || reportType == null || reportOutputPath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checks", options);
            System.exit(1);
        }

        RunContext runContext = null;
        try {
            runDirectory = FolderChecker.build().checkFolder(runDirectory);
            runContext = ProductionRunContextFactory.fromRunDirectory(runDirectory);
        } catch (IOException | HartwigException e) {
            LOGGER.info(e.getMessage());
            System.exit(1);
        }

        new HealthChecksApplication(runContext, reportType, reportOutputPath).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, "The path containing the data for a single run");
        options.addOption(REPORT_TYPE, true, "The type of report to be generated: 'json' or 'stdout'.");
        options.addOption(REPORT_OUTPUT_PATH, true, "The path where reports are written to.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private void run() {
        final HealthChecksFlyweight flyweight = HealthChecksFlyweight.getInstance();
        final Collection<HealthChecker> checkers = flyweight.getAllCheckers();

        final Observable<HealthChecker> checkerObservable = Observable.from(checkers).subscribeOn(Schedulers.io());

        BlockingObservable.from(checkerObservable).subscribe(createHealthCheckerAction(), createErrorHandler(),
                this::generateReport);
    }

    @NotNull
    private Action1<? super HealthChecker> createHealthCheckerAction() {
        return (Action1<HealthChecker>) healthChecker -> {
            final Report report = HealthCheckReportFactory.create(reportType);
            report.addResult(healthChecker.run(runContext));
        };
    }

    @NotNull
    private static Action1<? super Throwable> createErrorHandler() {
        return (Action1<Throwable>) throwable -> LOGGER.error(throwable.getMessage());
    }

    private void generateReport() {
        try {
            final Report report = HealthCheckReportFactory.create(reportType);

            final Optional<String> reportPath = report.generateReport(runContext, reportOutputPath);
            if (reportPath.isPresent()) {
                LOGGER.info(String.format("Report generated -> \n%s", reportPath.get()));
            }
        } catch (final GenerateReportException e) {
            LOGGER.log(Level.ERROR, e.getMessage());
        }
    }
}
