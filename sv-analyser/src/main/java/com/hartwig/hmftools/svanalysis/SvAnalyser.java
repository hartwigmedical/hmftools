package com.hartwig.hmftools.svanalysis;

import static com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator.PON_FILTER_PON;

import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.analysis.CNAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.FusionDisruptionAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.SvClusteringConfig;
import com.hartwig.hmftools.svanalysis.analysis.SvSampleAnalyser;
import com.hartwig.hmftools.svanalysis.annotators.ExtDataLinker;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser;

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


public class SvAnalyser {

    private static final Logger LOGGER = LogManager.getLogger(SvAnalyser.class);

    public static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String RUN_SVA = "run_sv_analysis";
    private static final String DATA_OUTPUT_PATH = "data_output_path";
    private static final String LOG_DEBUG = "log_debug";
    private static final String WRITE_FILTERED_SVS = "write_pon_filters";
    private static final String LOG_VCF_INSERTS = "log_vcf_inserts";
    private static final String LOG_VCF_MANTA_DATA = "log_vcf_manta_data";
    private static final String LINE_ELEMENT_FILE = "line_element_file";
    private static final String COPY_NUMBER_ANALYSIS = "run_cn_analysis";
    private static final String RUN_RESULTS_CHECKER = "run_results_checker";
    private static final String EXTERNAL_DATA_LINK_FILE = "ext_data_link_file";
    private static final String INCLUDE_NONE_SEGMENTS = "incl_none_segments";
    private static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";

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

        if(tumorSample == null || tumorSample.equals("*"))
            tumorSample = "";

        List<String> samplesList = Lists.newArrayList();

        if (tumorSample.isEmpty())
        {
            if(dbAccess != null)
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

        if(cmd.hasOption(COPY_NUMBER_ANALYSIS))
        {
            CNAnalyser cnAnalyser = new CNAnalyser(cmd.getOptionValue(DATA_OUTPUT_PATH), dbAccess);

            cnAnalyser.loadConfig(cmd, samplesList);
            cnAnalyser.analyseData();
            cnAnalyser.close();

            LOGGER.info("CN analysis complete");
            return;
        }

        if(cmd.hasOption(RUN_SVA))
        {
            ExtDataLinker extDataLinker = null;

            if (cmd.hasOption(EXTERNAL_DATA_LINK_FILE))
            {
                extDataLinker = new ExtDataLinker();
                extDataLinker.loadFile(cmd.getOptionValue(EXTERNAL_DATA_LINK_FILE));
            }

            SvClusteringConfig clusteringConfig = new SvClusteringConfig(cmd, tumorSample);
            SvSampleAnalyser sampleAnalyser = new SvSampleAnalyser(clusteringConfig);

            FusionDisruptionAnalyser fusionAnalyser = null;

            if(cmd.hasOption(GENE_TRANSCRIPTS_DIR))
            {
                fusionAnalyser = new FusionDisruptionAnalyser();
                fusionAnalyser.loadFusionReferenceData(cmd, cmd.getOptionValue(DATA_OUTPUT_PATH), samplesList.size() > 1);
            }

            int count = 0;
            for (final String sample : samplesList)
            {
                ++count;
                List<SvVarData> svVarData = queryStructuralVariantData(dbAccess, sample);

                if(svVarData.isEmpty())
                    continue;

                LOGGER.info("sample({}) processing {} SVs, totalProcessed({})", sample, svVarData.size(), count);

                if (cmd.hasOption(INCLUDE_NONE_SEGMENTS))
                {
                    CNAnalyser cnAnalyser = new CNAnalyser(cmd.getOptionValue(DATA_OUTPUT_PATH), dbAccess);

                    int varCount = svVarData.size();
                    List<StructuralVariantData> noneSegmentSVs = cnAnalyser.loadNoneSegments(sample, varCount + 1);

                    LOGGER.debug("sample({}) including {} none copy number segments", sample, noneSegmentSVs.size());

                    for (final StructuralVariantData svData : noneSegmentSVs)
                    {
                        SvVarData var = new SvVarData(svData);
                        var.setNoneSegment(true);
                        svVarData.add(var);
                    }
                }

                sampleAnalyser.loadFromDatabase(sample, svVarData);

                if (extDataLinker != null && extDataLinker.hasData())
                {
                    extDataLinker.setSVData(sample, svVarData);
                }
                else
                {
                    sampleAnalyser.analyse();

                    if(fusionAnalyser != null)
                    {
                        fusionAnalyser.loadSvGeneTranscriptData(sample, cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));
                        fusionAnalyser.findFusions(svVarData, sampleAnalyser.getClusters());
                    }
                }
            }

            sampleAnalyser.close();

            if(fusionAnalyser != null)
                fusionAnalyser.close();
        }

        if(cmd.hasOption(RUN_RESULTS_CHECKER))
        {
            ResultsChecker resultsChecker = new ResultsChecker();
            resultsChecker.setLogMismatches(cmd.hasOption(LOG_DEBUG));

            if(resultsChecker.loadConfig(cmd, samplesList, cmd.getOptionValue(DATA_OUTPUT_PATH)))
            {
                resultsChecker.setIdColumns(true);
                resultsChecker.addDefaultColumnsToCheck();

                if(resultsChecker.loadData())
                {
                    if (resultsChecker.runChecks())
                        LOGGER.info("results validation passed");
                    else
                        LOGGER.warn("results validation failed");
                }
            }
        }

        if(cmd.hasOption(WRITE_FILTERED_SVS) || cmd.hasOption(LOG_VCF_INSERTS))
        {
            LOGGER.info("reading VCF files including filtered SVs");

            FilteredSVWriter filteredSvWriter = new FilteredSVWriter(cmd.getOptionValue(VCF_FILE), cmd.getOptionValue(DATA_OUTPUT_PATH));

            filteredSvWriter.setLogInsSVs(cmd.hasOption(LOG_VCF_INSERTS));
            filteredSvWriter.setRunPONFilter(cmd.hasOption(WRITE_FILTERED_SVS));
            filteredSvWriter.setLogExtraMantaData(cmd.hasOption(LOG_VCF_MANTA_DATA));
            filteredSvWriter.processVcfFiles();

            LOGGER.info("reads complete");
            return;
        }

        LOGGER.info("run complete");
    }

    private static List<SvVarData> queryStructuralVariantData(@NotNull DatabaseAccess dbAccess, @NotNull String sampleId)
    {
        List<SvVarData> svVarDataItems = Lists.newArrayList();

        List<StructuralVariantData> svRecords = dbAccess.readStructuralVariantData(sampleId);

        for (final StructuralVariantData svRecord : svRecords) {

            if(svRecord.filter().equals(PON_FILTER_PON))
                continue;

            if(!svRecord.filter().equals("PASS"))
                continue;

            svVarDataItems.add(new SvVarData(svRecord));
        }

        return svVarDataItems;
    }

    private static List<String> getStructuralVariantSamplesList(@NotNull DatabaseAccess dbAccess)
    {
        return dbAccess.structuralVariantSampleList("");
    }

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(RUN_SVA, false, "Whether to run clustering logic");
        options.addOption(DATA_OUTPUT_PATH, true, "CSV output directory");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(WRITE_FILTERED_SVS, false, "Includes filtered SVs and writes all to file for PON creation");
        options.addOption(LOG_VCF_INSERTS, false, "Read INS from VCF files, write to CSV");
        options.addOption(LOG_VCF_MANTA_DATA, false, "Read extra manta data from VCF files, write to CSV");
        options.addOption(LINE_ELEMENT_FILE, true, "Line Elements file for SVs");
        options.addOption(COPY_NUMBER_ANALYSIS, false, "Run copy number analysis");
        options.addOption(RUN_RESULTS_CHECKER, false, "Check results vs validation file");
        options.addOption(INCLUDE_NONE_SEGMENTS, false, "Include copy number NONE segments in SV analysis");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Optinonal: file with sample gene transcript data");
        SvClusteringConfig.addCmdLineArgs(options);
        ResultsChecker.addCmdLineArgs(options);
        CNAnalyser.addCmdLineArgs(options);
        SvFusionAnalyser.addCmdLineArgs(options);

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
