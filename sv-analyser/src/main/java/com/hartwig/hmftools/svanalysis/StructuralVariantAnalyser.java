package com.hartwig.hmftools.svanalysis;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.analysis.StructuralVariantClustering;
import com.hartwig.hmftools.svanalysis.annotators.SvVCFAnnotator;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.analysis.SvClusteringConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;


public class StructuralVariantAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnalyser.class);

    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String CLUSTER_SVS = "cluster_svs";
    private static final String CLUSTER_BASE_DISTANCE = "cluster_bases";
    private static final String DATA_OUTPUT_PATH = "data_output_path";
    private static final String LOG_DEBUG = "log_debug";
    private static final String WRITE_FILTERED_SVS = "log_filters";
    private static final String SV_PON_FILE = "sv_pon_file";
    private static final String FRAGILE_SITE_FILE = "fragile_site_file";
    private static final String LINE_ELEMENT_FILE = "line_element_file";
    private static final String REANNOTATE_FROM_VCFS = "reannotate_vcfs";
    private static final String EXTERNAL_ANNOTATIONS = "external_annotations";


    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        if(cmd.hasOption(WRITE_FILTERED_SVS))
        {
            LOGGER.info("reading VCF files including filtered SVs");

            FilteredSVWriter filteredSvWriter = new FilteredSVWriter(cmd.getOptionValue(VCF_FILE), cmd.getOptionValue(DATA_OUTPUT_PATH));
            filteredSvWriter.processVcfFiles();

            LOGGER.info("reads complete");
            return;
        }

        if(cmd.hasOption(REANNOTATE_FROM_VCFS))
        {
            LOGGER.info("reading VCF files to re-annotate");

            // for now just re-read the VCFs and write out new annotations to file
            // may later on turn into update SQL once clustering does the same
            SvVCFAnnotator vcfAnnotator = new SvVCFAnnotator(cmd.getOptionValue(VCF_FILE), cmd.getOptionValue(DATA_OUTPUT_PATH));
            vcfAnnotator.processVcfFiles();
            return;
        }

        final DatabaseAccess dbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

        final String tumorSample = cmd.getOptionValue(SAMPLE);

        SvClusteringConfig clusteringConfig = new SvClusteringConfig();
        clusteringConfig.setOutputCsvPath(cmd.getOptionValue(DATA_OUTPUT_PATH));
        clusteringConfig.setBaseDistance(Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE, "0")));
        clusteringConfig.setUseCombinedOutputFile(tumorSample.equals("*"));
        clusteringConfig.setSvPONFile(cmd.getOptionValue(SV_PON_FILE, ""));
        clusteringConfig.setFragileSiteFile(cmd.getOptionValue(FRAGILE_SITE_FILE, ""));
        clusteringConfig.setLineElementFile(cmd.getOptionValue(LINE_ELEMENT_FILE, ""));
        clusteringConfig.setExternalAnnotationsFile(cmd.getOptionValue(EXTERNAL_ANNOTATIONS, ""));
        StructuralVariantClustering svClusterer = new StructuralVariantClustering(clusteringConfig);

        List<String> samplesList = Lists.newArrayList();

        if (tumorSample.isEmpty() || tumorSample.equals("*")) {
            samplesList = getStructuralVariantSamplesList(dbAccess);
        } else if (tumorSample.contains(",")) {
            String[] tumorList = tumorSample.split(",");
            samplesList = Arrays.stream(tumorList).collect(Collectors.toList());
        } else {
            samplesList.add(tumorSample);
        }

        int count = 0;
        for (final String sample : samplesList) {
            ++count;
            LOGGER.info("clustering for sample({}), total({})", sample, count);

            List<SvClusterData> svClusterData = queryStructuralVariantData(dbAccess, sample);

            svClusterer.loadFromDatabase(sample, svClusterData);

            //LOGGER.info("data loaded", sample, count);

            svClusterer.runClustering();
            //LOGGER.info("clustering complete", sample, count);

//                if(count > 10)
//                    break;
        }

        svClusterer.close();

        LOGGER.info("run complete");
    }

    @NotNull
    private static List<SvClusterData> queryStructuralVariantData(@NotNull DatabaseAccess dbAccess, @NotNull String sampleId) {
        List<SvClusterData> svClusterDataItems = Lists.newArrayList();

        List<StructuralVariantData> svRecords = dbAccess.readStructuralVariantData(sampleId);

        for (final StructuralVariantData svRecord : svRecords) {
            svClusterDataItems.add(new SvClusterData(svRecord));
        }

        return svClusterDataItems;
    }

    @NotNull
    private static List<String> getStructuralVariantSamplesList(@NotNull DatabaseAccess dbAccess) {
        return dbAccess.structuralVariantSampleList("");
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(CLUSTER_SVS, false, "Whether to run clustering logic");
        options.addOption(DATA_OUTPUT_PATH, true, "CSV output directory");
        options.addOption(CLUSTER_BASE_DISTANCE, true, "Clustering base distance, defaults to 1000");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(WRITE_FILTERED_SVS, false, "Includes filtered SVs and writes all to file for PON creation");
        options.addOption(SV_PON_FILE, true, "PON file for SVs");
        options.addOption(LINE_ELEMENT_FILE, true, "Line Elements file for SVs");
        options.addOption(FRAGILE_SITE_FILE, true, "Fragile Site file for SVs");
        options.addOption(REANNOTATE_FROM_VCFS, false, "Re-annotation with new VCF annotations");
        options.addOption(EXTERNAL_ANNOTATIONS, true, "Exteranl file with per-SV annotations");

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
