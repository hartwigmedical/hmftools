package com.hartwig.hmftools.svanalysis;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.NONE_SEGMENT_INFERRED;

import java.io.File;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.analysis.CNAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.FusionDisruptionAnalyser;
import com.hartwig.hmftools.svanalysis.stats.StatisticRoutines;
import com.hartwig.hmftools.svanalysis.analysis.SvaConfig;
import com.hartwig.hmftools.svanalysis.analysis.SvSampleAnalyser;
import com.hartwig.hmftools.svanalysis.annotators.DriverGeneAnnotator;
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
    private static final String LINE_ELEMENT_FILE = "line_element_file";
    private static final String DRIVERS_CHECK = "check_drivers";
    private static final String RUN_FUSIONS = "run_fusions";
    private static final String COPY_NUMBER_ANALYSIS = "run_cn_analysis";
    private static final String INCLUDE_NONE_SEGMENTS = "incl_none_segments";
    private static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    private static final String STATS_ROUTINES = "stats_routines";

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

        if(cmd.hasOption(STATS_ROUTINES))
        {
            StatisticRoutines statsRoutines = new StatisticRoutines();
            statsRoutines.loadConfig(cmd);
            statsRoutines.runStatistics();
            LOGGER.info("run complete");
            return;
        }

        final DatabaseAccess dbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

        String sampleId = cmd.getOptionValue(SAMPLE);

        if(sampleId == null || sampleId.equals("*"))
            sampleId = "";

        List<String> samplesList = Lists.newArrayList();

        if (sampleId.isEmpty())
        {
            samplesList = getStructuralVariantSamplesList(dbAccess);
        }
        else if (sampleId.contains(","))
        {
            String[] tumorList = sampleId.split(",");
            samplesList = Arrays.stream(tumorList).collect(Collectors.toList());
        }
        else
        {
            samplesList.add(sampleId);
        }

        String dataOutputDir = cmd.getOptionValue(DATA_OUTPUT_PATH, "");

        if(!dataOutputDir.endsWith(File.separator))
            dataOutputDir += File.separator;

        if(cmd.hasOption(COPY_NUMBER_ANALYSIS))
        {
            CNAnalyser cnAnalyser = new CNAnalyser(dataOutputDir, dbAccess);

            cnAnalyser.loadConfig(cmd, samplesList);
            cnAnalyser.findLOHEvents();
            cnAnalyser.close();

            LOGGER.info("CN analysis complete");
            return;
        }

        if(cmd.hasOption(RUN_SVA))
        {
            SvaConfig svaConfig = new SvaConfig(cmd, sampleId);
            SvSampleAnalyser sampleAnalyser = new SvSampleAnalyser(svaConfig);

            DriverGeneAnnotator driverGeneAnnotator = null;
            boolean checkDrivers = cmd.hasOption(DRIVERS_CHECK);

            FusionDisruptionAnalyser fusionAnalyser = null;
            boolean runFusions = cmd.hasOption(RUN_FUSIONS);

            if((runFusions || checkDrivers) && cmd.hasOption(GENE_TRANSCRIPTS_DIR))
            {
                fusionAnalyser = new FusionDisruptionAnalyser();
                fusionAnalyser.loadFusionReferenceData(cmd, dataOutputDir, cmd.getOptionValue(GENE_TRANSCRIPTS_DIR,""));

                if(!fusionAnalyser.getRnaSampleIds().isEmpty())
                {
                    samplesList.clear();
                    samplesList.addAll(fusionAnalyser.getRnaSampleIds());

                    LOGGER.info("running {} sample based on RNA fusion input", samplesList.size());
                }

                // fusionAnalyser.writeGeneProbabilityData();
            }

            if(checkDrivers)
            {
                driverGeneAnnotator = new DriverGeneAnnotator(dbAccess, fusionAnalyser.getGeneTranscriptCollection(), dataOutputDir);
                driverGeneAnnotator.loadConfig(cmd);
            }

            CNAnalyser cnAnalyser = new CNAnalyser(dataOutputDir, dbAccess);
            boolean createNoneSvsFromCNData = cmd.hasOption(INCLUDE_NONE_SEGMENTS);

            if(!svaConfig.LOHDataFile.isEmpty())
            {
                cnAnalyser.loadLOHFromCSV(svaConfig.LOHDataFile, "");
                sampleAnalyser.setSampleLohData(cnAnalyser.getSampleLohData());

                if(driverGeneAnnotator != null)
                    driverGeneAnnotator.setLohData(cnAnalyser.getSampleLohData());
            }

            int count = 0;
            for (final String sample : samplesList)
            {
                ++count;
                List<SvVarData> svVarData = queryStructuralVariantData(dbAccess, sample, !createNoneSvsFromCNData);

                if(svVarData.isEmpty())
                {
                    LOGGER.debug("sample({}) has no SVs, totalProcessed({})", sample, count);
                    continue;
                }

                LOGGER.info("sample({}) processing {} SVs, totalProcessed({})", sample, svVarData.size(), count);

                List<SegmentSupport> cnSegmentTypes = Lists.newArrayList();

                cnSegmentTypes.add(SegmentSupport.CENTROMERE);
                cnSegmentTypes.add(SegmentSupport.TELOMERE);

                if (createNoneSvsFromCNData)
                    cnSegmentTypes.add(SegmentSupport.NONE);

                final List<PurpleCopyNumber> cnRecords = dbAccess.readCopyNumberSegmentsByType(sample, cnSegmentTypes);

                final Map<String, double[]> chrCopyNumberMap = cnAnalyser.createChrCopyNumberMap(cnRecords);

                if (createNoneSvsFromCNData)
                {
                    List<StructuralVariantData> noneSegmentSVs = cnAnalyser.createNoneSegments(cnRecords);

                    LOGGER.debug("sample({}) including {} none copy number segments", sample, noneSegmentSVs.size());

                    for (final StructuralVariantData svData : noneSegmentSVs)
                    {
                        svVarData.add(new SvVarData(svData));
                    }
                }

                sampleAnalyser.loadFromDatabase(sample, svVarData);

                sampleAnalyser.setChrCopyNumberMap(chrCopyNumberMap);

                sampleAnalyser.analyse();

                if(runFusions || checkDrivers)
                {
                    fusionAnalyser.setSvGeneData(sample, svVarData, runFusions);
                }

                if(checkDrivers)
                {
                    final PurityContext purityContext = dbAccess.readPurityContext(sample);

                    if(purityContext != null)
                        driverGeneAnnotator.setSamplePloidy(purityContext.bestFit().ploidy());

                    driverGeneAnnotator.setChrCopyNumberMap(chrCopyNumberMap);
                    driverGeneAnnotator.annotateSVs(sample, sampleAnalyser.getClusters(), sampleAnalyser.getChrBreakendMap());
                }

                sampleAnalyser.writeOutput();

                if(runFusions)
                {
                    fusionAnalyser.findFusions(svVarData, sampleAnalyser.getClusters());
                }

                if(svaConfig.MaxSamples > 0 && count >= svaConfig.MaxSamples)
                {
                    LOGGER.info("exiting after max sample count {} reached", count);
                    break;
                }
            }

            sampleAnalyser.close();

            if(fusionAnalyser != null)
                fusionAnalyser.close();

            if(driverGeneAnnotator != null)
                driverGeneAnnotator.close();

            LOGGER.info("SV analysis complete");
        }

        LOGGER.info("run complete");
    }

    private static List<SvVarData> queryStructuralVariantData(final DatabaseAccess dbAccess,
            final String sampleId, final boolean includeNoneSegments)
    {
        List<SvVarData> svVarDataItems = Lists.newArrayList();

        List<StructuralVariantData> svRecords = dbAccess.readStructuralVariantData(sampleId);

        for (final StructuralVariantData svRecord : svRecords) {

            if(svRecord.filter().equals(PON_FILTER_PON))
                continue;
            else if(svRecord.filter().equals(NONE_SEGMENT_INFERRED) && !includeNoneSegments)
                continue;

            // all others (currently PASS or blank) are accepted
            svVarDataItems.add(new SvVarData(svRecord));
        }

        return svVarDataItems;
    }

    private static List<String> getStructuralVariantSamplesList(@NotNull DatabaseAccess dbAccess)
    {
        final List<String> sampleIds = dbAccess.getSamplesPassingQC(MIN_SAMPLE_PURITY);

        if(!sampleIds.isEmpty())
            return sampleIds;

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
        options.addOption(DRIVERS_CHECK, false, "Check SVs against drivers catalog");
        options.addOption(RUN_FUSIONS, false, "Run fusion detection");
        options.addOption(COPY_NUMBER_ANALYSIS, false, "Run copy number analysis");
        options.addOption(INCLUDE_NONE_SEGMENTS, false, "Include copy number NONE segments in SV analysis");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Optional: file with sample gene transcript data");
        options.addOption(STATS_ROUTINES, false, "Optional: calc stats routines");
        SvaConfig.addCmdLineArgs(options);
        CNAnalyser.addCmdLineArgs(options);
        SvFusionAnalyser.addCmdLineArgs(options);
        StatisticRoutines.addCmdLineArgs(options);
        DriverGeneAnnotator.addCmdLineArgs(options);

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
