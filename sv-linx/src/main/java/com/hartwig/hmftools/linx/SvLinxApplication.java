package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;
import static com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser.setSvGeneData;
import static com.hartwig.hmftools.linx.types.SvaConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.linx.types.SvaConfig.LOG_DEBUG;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.types.SvaConfig;
import com.hartwig.hmftools.linx.analysis.SvSampleAnalyser;
import com.hartwig.hmftools.linx.drivers.DriverGeneAnnotator;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.fusion.FusionFinder;

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


public class SvLinxApplication
{

    private static final String DRIVERS_CHECK = "check_drivers";
    private static final String CHECK_FUSIONS = "check_fusions";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String FILTER_QC_PASS = "filter_qc_pass";

    private static final Logger LOGGER = LogManager.getLogger(SvLinxApplication.class);

    public static void main(@NotNull final String[] args) throws ParseException, SQLException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        SvaConfig svaConfig = new SvaConfig(cmd);

        final DatabaseAccess dbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

        List<String> samplesList = svaConfig.getSampleIds();

        if(dbAccess == null && (svaConfig.hasMultipleSamples() || samplesList.isEmpty()))
        {
            LOGGER.warn("batch mode requires a DB connection");
            return;
        }

        boolean sampleDataFromFile = (dbAccess == null);

        if (samplesList.isEmpty())
        {
            boolean filterQCPassOnly = cmd.hasOption(FILTER_QC_PASS);
            samplesList = getStructuralVariantSamplesList(dbAccess, filterQCPassOnly);

            LOGGER.info("retrieved {} samples {}", samplesList.size(), filterQCPassOnly ? "QC-pass-filtered" : "");

            svaConfig.setSampleIds(samplesList);
        }

        CnDataLoader cnDataLoader = new CnDataLoader(svaConfig.PurpleDataPath, dbAccess);

        SvSampleAnalyser sampleAnalyser = new SvSampleAnalyser(svaConfig);

        sampleAnalyser.setCnDataLoader(cnDataLoader);

        DriverGeneAnnotator driverGeneAnnotator = null;
        boolean checkDrivers = cmd.hasOption(DRIVERS_CHECK);

        FusionDisruptionAnalyser fusionAnalyser = null;
        boolean checkFusions = cmd.hasOption(CHECK_FUSIONS);

        boolean selectiveGeneLoading = (samplesList.size() == 1) && !checkDrivers;

        SvGeneTranscriptCollection ensemblDataCache = null;

        if(cmd.hasOption(GENE_TRANSCRIPTS_DIR))
        {
            ensemblDataCache = new SvGeneTranscriptCollection();
            ensemblDataCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));

            if(!ensemblDataCache.loadEnsemblData(selectiveGeneLoading))
            {
                LOGGER.error("Ensembl data cache load failed, exiting");
                return;
            }

            sampleAnalyser.setGeneCollection(ensemblDataCache);
            sampleAnalyser.getVisWriter().setGeneDataCollection(ensemblDataCache);

            if(checkFusions)
            {
                fusionAnalyser = new FusionDisruptionAnalyser();

                fusionAnalyser.initialise(cmd, svaConfig.OutputDataPath, svaConfig, ensemblDataCache);
                fusionAnalyser.setVisWriter(sampleAnalyser.getVisWriter());

                if(fusionAnalyser.hasRnaSampleData() && samplesList.size() > 1)
                {
                    samplesList.clear();
                    samplesList.addAll(fusionAnalyser.getRnaSampleIds());

                    LOGGER.info("running {} sample based on RNA fusion input", samplesList.size());
                }
            }

            if(checkDrivers)
            {
                driverGeneAnnotator = new DriverGeneAnnotator(dbAccess, ensemblDataCache, svaConfig.OutputDataPath);
                driverGeneAnnotator.loadConfig(cmd);
                driverGeneAnnotator.setVisWriter(sampleAnalyser.getVisWriter());
            }
        }

        if(driverGeneAnnotator != null)
        {
            driverGeneAnnotator.setCopyNumberData(
                    cnDataLoader.getChrCopyNumberMap(),
                    cnDataLoader.getLohData(),
                    cnDataLoader.getHomLossData());
        }

        PerformanceCounter prefCounter = new PerformanceCounter("SVA Total");

        int count = 0;
        for (final String sampleId : samplesList)
        {
            ++count;

            prefCounter.start();

            List<StructuralVariantData> svRecords = sampleDataFromFile ?
                    loadSampleSvData(svaConfig.SvDataPath, sampleId) : dbAccess.readStructuralVariantData(sampleId);

            List<SvVarData> svVarData = createSvData(svRecords);

            if(svVarData.isEmpty())
            {
                LOGGER.debug("sample({}) has no SVs, totalProcessed({})", sampleId, count);
                continue;
            }

            if(svaConfig.hasMultipleSamples())
            {
                LOGGER.info("sample({}) processing {} SVs, samplesComplete({})", sampleId, svVarData.size(), count);
            }

            cnDataLoader.loadSampleData(sampleId, svRecords);

            sampleAnalyser.setSampleSVs(sampleId, svVarData);

            sampleAnalyser.analyse();

            if(!sampleAnalyser.inValidState())
            {
                LOGGER.info("exiting after sample({}), in invalid state", sampleId);
                break;
            }

            if(ensemblDataCache != null)
            {
                // when matching RNA, allow all transcripts regardless of their viability for fusions
                boolean keepInvalidTranscripts = fusionAnalyser != null && fusionAnalyser.hasRnaSampleData();
                setSvGeneData(svVarData, ensemblDataCache, true, selectiveGeneLoading, !keepInvalidTranscripts);

                sampleAnalyser.annotateWithGeneData(ensemblDataCache);
            }

            if(checkDrivers)
            {
                driverGeneAnnotator.annotateSVs(sampleId, sampleAnalyser.getClusters(), sampleAnalyser.getChrBreakendMap());
            }

            if(checkFusions)
            {
                fusionAnalyser.run(sampleId, svVarData, dbAccess,
                        sampleAnalyser.getClusters(), sampleAnalyser.getChrBreakendMap());
            }

            sampleAnalyser.writeOutput(dbAccess);

            prefCounter.stop();

            if(svaConfig.MaxSamples > 0 && count >= svaConfig.MaxSamples)
            {
                LOGGER.info("exiting after max sample count {} reached", count);
                break;
            }
        }

        if(LOGGER.isDebugEnabled() || svaConfig.hasMultipleSamples())
        {
            prefCounter.logStats();
        }

        sampleAnalyser.close();

        if(fusionAnalyser != null)
            fusionAnalyser.close();

        if(driverGeneAnnotator != null)
            driverGeneAnnotator.close();

        LOGGER.info("SV analysis complete");
    }

    private static List<StructuralVariantData> loadSampleSvData(final String samplePath, final String sampleId)
    {
        try
        {
            final String svDataFile = StructuralVariantFile.generateFilename(samplePath, sampleId);
            return StructuralVariantFile.read(svDataFile);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load SV data: {}", e.toString());
            return Lists.newArrayList();
        }
    }

    private static List<SvVarData> createSvData(List<StructuralVariantData> svRecords)
    {
        List<SvVarData> svVarDataItems = Lists.newArrayList();

        for (final StructuralVariantData svRecord : svRecords) {

            if(svRecord.filter().equals(PON_FILTER_PON))
                continue;

            // all others (currently PASS or blank) are accepted
            svVarDataItems.add(new SvVarData(svRecord));
        }

        return svVarDataItems;
    }


    private static List<String> getStructuralVariantSamplesList(@NotNull DatabaseAccess dbAccess, boolean filterQCPassOnly)
    {
        final List<String> sampleIds = filterQCPassOnly ? dbAccess.getSamplesPassingQC(MIN_SAMPLE_PURITY) : dbAccess.getSampleIds();

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
        options.addOption(DRIVERS_CHECK, false, "Check SVs against drivers catalog");
        options.addOption(CHECK_FUSIONS, false, "Run fusion detection");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Optional: Ensembl data cache directory");
        options.addOption(FILTER_QC_PASS, false, "Optional: If present will filter out QC-fail sample");

        // allow sub-components to add their specific config
        SvaConfig.addCmdLineArgs(options);
        FusionFinder.addCmdLineArgs(options);
        DriverGeneAnnotator.addCmdLineArgs(options);
        FusionDisruptionAnalyser.addCmdLineArgs(options);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException
    {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
