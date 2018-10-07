package com.hartwig.hmftools.svanalysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.*;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.analysis.CNAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.SvClusteringConfig;
import com.hartwig.hmftools.svanalysis.analysis.SvSampleAnalyser;
import com.hartwig.hmftools.svanalysis.annotators.ExtDataLinker;
import com.hartwig.hmftools.svanalysis.svgraph.BreakpointGraph;
import com.hartwig.hmftools.svanalysis.svgraph.SimpleSimplificationStrategy;
import com.hartwig.hmftools.svanalysis.svgraph.Simplification;
import com.hartwig.hmftools.svanalysis.svgraph.SimplificationFile;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator.PON_FILTER_PON;


public class GraphAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(GraphAnalyser.class);

    private static final String LOG_DEBUG = "log_debug";
    private static final String SAMPLE = "sample";
    private static final String OUTPUT_DIRECTORY = "output_dir";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, SQLException, IOException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        final DatabaseAccess dbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

        String tumorSample = cmd.getOptionValue(SAMPLE);

        if(tumorSample == null || tumorSample.equals("*")) {
            tumorSample = "";
        }

        SvClusteringConfig clusteringConfig = new SvClusteringConfig(cmd, tumorSample);

        List<String> samplesList = Lists.newArrayList();

        if (tumorSample.isEmpty())
        {
            samplesList = getStructuralVariantSamplesList(dbAccess);
        }
        else if (tumorSample.contains(","))
        {
            String[] tumorList = tumorSample.split(",");
            samplesList = Arrays.stream(tumorList).collect(Collectors.toList());
        }
        else
        {
            samplesList.add(tumorSample);
        }

        for (final String sample : samplesList) {
            List<EnrichedStructuralVariant> svRecords = dbAccess.readStructuralVariants(sample);
            svRecords = hackFixPurpleOffByOne(svRecords);
            List<PurpleCopyNumber> cnRecords = dbAccess.readCopynumbers(sample);

            LOGGER.info("sample({}) processing {} SVs, {} segments", sample, svRecords.size(), cnRecords.size());

            BreakpointGraph graph = new BreakpointGraph(cnRecords, svRecords);
            List<Simplification> simplifications = graph.simplify(new SimpleSimplificationStrategy());

            String cnFile = PurpleCopyNumberFile.generateFilename(cmd.getOptionValue(OUTPUT_DIRECTORY), String.format("%s_cn_reduced", sample));
            String remainingFile = EnrichedStructuralVariantFile.generateFilename(cmd.getOptionValue(OUTPUT_DIRECTORY), String.format("%s_sv_remaining.csv", sample));
            String reducedFile = EnrichedStructuralVariantFile.generateFilename(cmd.getOptionValue(OUTPUT_DIRECTORY), String.format("%s_sv_simplified.csv", sample));
            String simplificationFile = SimplificationFile.generateFilename(cmd.getOptionValue(OUTPUT_DIRECTORY), String.format("%s", sample));

            PurpleCopyNumberFile.write(cnFile, graph.getAllSegments().stream().map(s -> s.cn()).collect(Collectors.toList()));
            EnrichedStructuralVariantFile.write(remainingFile, graph.getAllStructuralVariants());
            EnrichedStructuralVariantFile.write(reducedFile, simplifications.stream().flatMap(s -> s.variants().stream()).collect(Collectors.toList()));
            SimplificationFile.write(simplificationFile, simplifications);
        }
        LOGGER.info("run complete");
    }

    private static List<EnrichedStructuralVariant> hackFixPurpleOffByOne(List<EnrichedStructuralVariant> svRecords) {
        return svRecords.stream().map(sv ->
            ImmutableEnrichedStructuralVariant.builder().from(sv)
                    .start(hackFixPurpleOffByOne(sv.start()))
                    .end(hackFixPurpleOffByOne(sv.end()))
                    .build()
            ).collect(Collectors.toList());
    }

    private static EnrichedStructuralVariantLeg hackFixPurpleOffByOne(EnrichedStructuralVariantLeg leg) {
        if (leg == null || leg.orientation() == -1) return leg;
        return ImmutableEnrichedStructuralVariantLeg.builder().from(leg)
                .position(leg.position() - 1)
                .build();
    }

    private static List<String> getStructuralVariantSamplesList(@NotNull DatabaseAccess dbAccess)
    {
        return dbAccess.structuralVariantSampleList("");
    }

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(OUTPUT_DIRECTORY, true, "output directory containing reduced breakpoint graphs");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        SvClusteringConfig.addCmdLineArgs(options);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
