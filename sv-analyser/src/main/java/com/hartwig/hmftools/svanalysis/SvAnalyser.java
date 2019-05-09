package com.hartwig.hmftools.svanalysis;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;
import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.CN_ANALYSIS_ONLY;
import static com.hartwig.hmftools.svanalysis.analysis.FusionDisruptionAnalyser.setSvGeneData;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.NONE_SEGMENT_INFERRED;
import static com.hartwig.hmftools.svanalysis.types.SvaConfig.DATA_OUTPUT_PATH;
import static com.hartwig.hmftools.svanalysis.types.SvaConfig.LOG_DEBUG;

import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.analysis.CNAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.FusionDisruptionAnalyser;
import com.hartwig.hmftools.svanalysis.simulation.SvSimulator;
import com.hartwig.hmftools.svanalysis.stats.MultipleBiopsyAnalyser;
import com.hartwig.hmftools.svanalysis.stats.StatisticRoutines;
import com.hartwig.hmftools.svanalysis.types.SvaConfig;
import com.hartwig.hmftools.svanalysis.analysis.SvSampleAnalyser;
import com.hartwig.hmftools.svanalysis.annotators.DriverGeneAnnotator;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;
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


    private static final String RUN_SVA = "run_sv_analysis";
    private static final String DRIVERS_CHECK = "check_drivers";
    private static final String CHECK_FUSIONS = "check_fusions";
    private static final String INCLUDE_NONE_SEGMENTS = "incl_none_segments";
    private static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    private static final String STATS_ROUTINES = "stats_routines";
    private static final String MULT_BIOPSY_ANALYSIS = "mult_biopsy_analysis";
    private static final String SIM_ROUTINES = "sim_routines";

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

        SvaConfig svaConfig = new SvaConfig(cmd);

        if(cmd.hasOption(STATS_ROUTINES))
        {
            StatisticRoutines statsRoutines = new StatisticRoutines();
            statsRoutines.loadConfig(cmd);
            statsRoutines.runStatistics();
            LOGGER.info("run complete");
            return;
        }

        if(cmd.hasOption(SIM_ROUTINES))
        {
            SvSimulator simulator = new SvSimulator();
            simulator.loadConfig(cmd, svaConfig.OutputCsvPath);
            simulator.run();
            return;
        }

        if(cmd.hasOption(MULT_BIOPSY_ANALYSIS))
        {
            MultipleBiopsyAnalyser mbAnalyser = new MultipleBiopsyAnalyser();

            if(mbAnalyser.loadData(cmd, svaConfig.OutputCsvPath))
            {
                mbAnalyser.runAnalysis();
            }

            return;
        }

        final DatabaseAccess dbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

        List<String> samplesList = Lists.newArrayList();

        if (svaConfig.SampleId.isEmpty())
        {
            samplesList = getStructuralVariantSamplesList(dbAccess);
        }
        else if (svaConfig.SampleId.contains(","))
        {
            String[] tumorList = svaConfig.SampleId.split(",");
            samplesList = Arrays.stream(tumorList).collect(Collectors.toList());
        }
        else
        {
            samplesList.add(svaConfig.SampleId);
        }

        CNAnalyser cnAnalyser = new CNAnalyser(svaConfig.OutputCsvPath, dbAccess);
        cnAnalyser.loadConfig(cmd, samplesList);

        if(cmd.hasOption(CN_ANALYSIS_ONLY))
        {
            // run CN analysis, which will write a bunch of cohort-wide sample data, then exit
            cnAnalyser.runSamplesList();
            cnAnalyser.close();

            LOGGER.info("CN analysis complete");
            return;
        }

        if(cmd.hasOption(RUN_SVA))
        {
            SvSampleAnalyser sampleAnalyser = new SvSampleAnalyser(svaConfig);

            sampleAnalyser.setCopyNumberAnalyser(cnAnalyser);

            DriverGeneAnnotator driverGeneAnnotator = null;
            boolean checkDrivers = cmd.hasOption(DRIVERS_CHECK);

            FusionDisruptionAnalyser fusionAnalyser = null;
            boolean checkFusions = cmd.hasOption(CHECK_FUSIONS);

            boolean isSingleSample = samplesList.size() == 1;

            SvGeneTranscriptCollection ensemblDataCache = null;

            if(cmd.hasOption(GENE_TRANSCRIPTS_DIR))
            {
                ensemblDataCache = new SvGeneTranscriptCollection();
                ensemblDataCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));

                if(!ensemblDataCache.loadEnsemblData(isSingleSample))
                {
                    LOGGER.error("Ensembl data cache load failed, exiting");
                    return;
                }

                sampleAnalyser.setGeneCollection(ensemblDataCache);
                sampleAnalyser.getVisWriter().setGeneDataCollection(ensemblDataCache);

                if(checkFusions)
                {
                    fusionAnalyser = new FusionDisruptionAnalyser();

                    fusionAnalyser.initialise(cmd, svaConfig.OutputCsvPath, ensemblDataCache);
                    fusionAnalyser.setVisWriter(sampleAnalyser.getVisWriter());

                    if(!fusionAnalyser.getRnaSampleIds().isEmpty() && samplesList.size() > 1)
                    {
                        samplesList.clear();
                        samplesList.addAll(fusionAnalyser.getRnaSampleIds());

                        LOGGER.info("running {} sample based on RNA fusion input", samplesList.size());
                    }

                    // fusionAnalyser.writeGeneProbabilityData();

                }

                if(checkDrivers)
                {
                    driverGeneAnnotator = new DriverGeneAnnotator(dbAccess, ensemblDataCache, svaConfig.OutputCsvPath);
                    driverGeneAnnotator.loadConfig(cmd);
                    driverGeneAnnotator.setVisWriter(sampleAnalyser.getVisWriter());
                }
            }

            boolean createNoneSvsFromCNData = cmd.hasOption(INCLUDE_NONE_SEGMENTS);

            if(driverGeneAnnotator != null)
            {
                driverGeneAnnotator.setLohData(cnAnalyser.getSampleLohData());
                driverGeneAnnotator.setChrCopyNumberMap(cnAnalyser.getChrCopyNumberMap());
            }

            PerformanceCounter prefCounter = new PerformanceCounter("SVA Total");

            int count = 0;
            for (final String sampleId : samplesList)
            {
                ++count;

                prefCounter.start();

                List<StructuralVariantData> svRecords = dbAccess.readStructuralVariantData(sampleId);

                List<SvVarData> svVarData = createSvData(svRecords, !createNoneSvsFromCNData);

                if(svVarData.isEmpty())
                {
                    LOGGER.debug("sample({}) has no SVs, totalProcessed({})", sampleId, count);
                    continue;
                }

                LOGGER.info("sample({}) processing {} SVs, samplesComplete({})", sampleId, svVarData.size(), count);

                cnAnalyser.loadSampleData(sampleId, svRecords, createNoneSvsFromCNData);

                if (createNoneSvsFromCNData)
                {
                    List<StructuralVariantData> noneSegmentSVs = cnAnalyser.createNoneSegments();

                    LOGGER.debug("sample({}) including {} none copy number segments", sampleId, noneSegmentSVs.size());

                    for (final StructuralVariantData svData : noneSegmentSVs)
                    {
                        svVarData.add(new SvVarData(svData));
                    }
                }

                sampleAnalyser.loadFromDatabase(sampleId, svVarData);

                sampleAnalyser.analyse();

                if(ensemblDataCache != null)
                {
                    // when matching RNA, allow all transcripts regardless of their viability for fusions
                    setSvGeneData(svVarData, ensemblDataCache, checkFusions, isSingleSample, fusionAnalyser.getRnaSampleIds().isEmpty());
                }

                if(checkDrivers)
                {
                    final PurityContext purityContext = dbAccess.readPurityContext(sampleId);

                    if(purityContext != null)
                        driverGeneAnnotator.setSamplePloidy(purityContext.bestFit().ploidy());

                    driverGeneAnnotator.annotateSVs(sampleId, sampleAnalyser.getClusters(), sampleAnalyser.getChrBreakendMap());
                }

                if(checkFusions)
                {
                    fusionAnalyser.run(sampleId, svVarData, sampleAnalyser.getClusters(), sampleAnalyser.getChrBreakendMap());
                }

                sampleAnalyser.writeOutput();

                prefCounter.stop();

                if(svaConfig.MaxSamples > 0 && count >= svaConfig.MaxSamples)
                {
                    LOGGER.info("exiting after max sample count {} reached", count);
                    break;
                }
            }

            prefCounter.logStats();

            sampleAnalyser.close();

            if(fusionAnalyser != null)
                fusionAnalyser.close();

            if(driverGeneAnnotator != null)
                driverGeneAnnotator.close();

            LOGGER.info("SV analysis complete");
        }

        LOGGER.info("run complete");
    }

    private static List<SvVarData> createSvData(List<StructuralVariantData> svRecords, boolean includeNoneSegments)
    {
        List<SvVarData> svVarDataItems = Lists.newArrayList();

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
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(RUN_SVA, false, "Whether to run clustering logic");
        options.addOption(DRIVERS_CHECK, false, "Check SVs against drivers catalog");
        options.addOption(CHECK_FUSIONS, false, "Run fusion detection");
        options.addOption(INCLUDE_NONE_SEGMENTS, false, "Include copy number NONE segments in SV analysis");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Optional: file with sample gene transcript data");
        options.addOption(STATS_ROUTINES, false, "Optional: calc stats routines");
        options.addOption(SIM_ROUTINES, false, "Optional: simulation routines");
        options.addOption(MULT_BIOPSY_ANALYSIS, false, "Optional: run multiple biopsy analysis");

        // allow sub-components to add their specific config
        SvSimulator.addCmdLineArgs(options);
        SvaConfig.addCmdLineArgs(options);
        CNAnalyser.addCmdLineArgs(options);
        SvFusionAnalyser.addCmdLineArgs(options);
        StatisticRoutines.addCmdLineArgs(options);
        DriverGeneAnnotator.addCmdLineArgs(options);
        FusionDisruptionAnalyser.addCmdLineArgs(options);
        MultipleBiopsyAnalyser.addCmdLineArgs(options);

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
