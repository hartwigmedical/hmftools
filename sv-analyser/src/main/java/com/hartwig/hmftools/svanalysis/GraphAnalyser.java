package com.hartwig.hmftools.svanalysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.analysis.CNAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.SvClusteringConfig;
import com.hartwig.hmftools.svanalysis.analysis.SvSampleAnalyser;
import com.hartwig.hmftools.svanalysis.annotators.ExtDataLinker;
import com.hartwig.hmftools.svanalysis.svgraph.BreakpointGraph;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import org.apache.commons.cli.*;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

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
    private static final String OUTPUT_OPERATION = "output_ops";
    private static final String OUTPUT_COPY_NUMBER = "output_cn";
    private static final String OUTPUT_STRUCTURAL_VARIANTS = "output_sv";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, SQLException
    {
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
        SvSampleAnalyser sampleAnalyser = new SvSampleAnalyser(clusteringConfig);

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
            List<PurpleCopyNumber> cnRecords = dbAccess.readCopynumbers(sample);

            LOGGER.info("sample({}) processing {} SVs, {} segments", sample, svRecords.size(), cnRecords.size());

            BreakpointGraph graph = new BreakpointGraph(cnRecords, svRecords);
            List<StructuralVariant> list = graph.simplify();

            throw new RuntimeException("NYI: write simplifications & reduced genome");
        }
        LOGGER.info("run complete");
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
        options.addOption(OUTPUT_OPERATION, true, "output file containing operations performed");
        options.addOption(OUTPUT_COPY_NUMBER, true, "Final copy number segmentation");
        options.addOption(OUTPUT_STRUCTURAL_VARIANTS, true, "Remaining unsimplified SVs");
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
