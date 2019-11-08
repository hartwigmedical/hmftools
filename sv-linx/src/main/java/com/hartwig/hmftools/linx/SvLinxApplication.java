package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.linx.LinxConfig.DB_PASS;
import static com.hartwig.hmftools.linx.LinxConfig.DB_URL;
import static com.hartwig.hmftools.linx.LinxConfig.DB_USER;
import static com.hartwig.hmftools.linx.LinxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.linx.LinxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.linx.LinxConfig.LOG_VERBOSE;
import static com.hartwig.hmftools.linx.LinxConfig.REF_GENOME_FILE;
import static com.hartwig.hmftools.linx.LinxConfig.databaseAccess;
import static com.hartwig.hmftools.linx.SvDataLoader.VCF_FILE;
import static com.hartwig.hmftools.linx.SvDataLoader.loadSvDataFromSvFile;
import static com.hartwig.hmftools.linx.SvDataLoader.loadSvDataFromVcf;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.version.VersionInfo;
import com.hartwig.hmftools.linx.analysis.SvSampleAnalyser;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.drivers.DriverGeneAnnotator;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.fusion.FusionFinder;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

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
    private static final String FILTER_QC_PASS = "filter_qc_pass";

    private static final Logger LOGGER = LogManager.getLogger(SvLinxApplication.class);

    public static void main(@NotNull final String[] args) throws ParseException, SQLException
    {
        final VersionInfo version = new VersionInfo("linx.version");
        LOGGER.info("LINX version: {}", version.version());

        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_VERBOSE))
        {
            Configurator.setRootLevel(Level.TRACE);
        }
        else if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        LinxConfig config = new LinxConfig(cmd);

        if(!config.hasValidPaths())
        {
            LOGGER.warn("invalid config paths");
            return;
        }

        final DatabaseAccess dbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

        boolean sampleDataFromFile = !config.PurpleDataPath.isEmpty();

        List<String> samplesList = config.getSampleIds();

        if(dbAccess == null)
        {
            if((config.hasMultipleSamples() || samplesList.isEmpty()))
            {
                LOGGER.warn("batch mode requires a DB connection");
                return;
            }

            if(!sampleDataFromFile)
            {
                LOGGER.warn("no DB connection and no purple data directory specified");
                return;
            }
        }


        if (samplesList.isEmpty())
        {
            boolean filterQCPassOnly = cmd.hasOption(FILTER_QC_PASS);
            samplesList = getStructuralVariantSamplesList(dbAccess, filterQCPassOnly);

            LOGGER.info("retrieved {} samples {}", samplesList.size(), filterQCPassOnly ? "QC-pass-filtered" : "");

            config.setSampleIds(samplesList);
        }

        CnDataLoader cnDataLoader = new CnDataLoader(config.PurpleDataPath, dbAccess);

        SvSampleAnalyser sampleAnalyser = new SvSampleAnalyser(config);

        sampleAnalyser.setCnDataLoader(cnDataLoader);

        DriverGeneAnnotator driverGeneAnnotator = null;
        boolean checkDrivers = cmd.hasOption(DRIVERS_CHECK);

        FusionDisruptionAnalyser fusionAnalyser = null;
        boolean checkFusions = cmd.hasOption(CHECK_FUSIONS);

        boolean selectiveGeneLoading = (samplesList.size() == 1) && !checkDrivers;
        boolean applyPromotorDistance = checkFusions;
        boolean purgeInvalidTranscripts = true;

        SvGeneTranscriptCollection ensemblDataCache = null;

        if(cmd.hasOption(GENE_TRANSCRIPTS_DIR))
        {
            ensemblDataCache = new SvGeneTranscriptCollection();
            ensemblDataCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));

            if(checkDrivers && !checkFusions)
                ensemblDataCache.setRequireCodingInfo(false);

            if(!ensemblDataCache.loadEnsemblData(selectiveGeneLoading))
            {
                LOGGER.error("Ensembl data cache load failed, exiting");
                return;
            }

            sampleAnalyser.setGeneCollection(ensemblDataCache);
            sampleAnalyser.getVisWriter().setGeneDataCollection(ensemblDataCache);

            // always initialise since is used for transcript evaluation
            fusionAnalyser = new FusionDisruptionAnalyser(
                    cmd, config, ensemblDataCache, sampleAnalyser.getVisWriter());

            if(checkFusions)
            {
                purgeInvalidTranscripts = !fusionAnalyser.hasRnaSampleData();

                if(fusionAnalyser.hasRnaSampleData() && samplesList.size() > 1)
                {
                    samplesList.clear();
                    samplesList.addAll(fusionAnalyser.getRnaSampleIds());

                    LOGGER.info("running {} sample based on RNA fusion input", samplesList.size());
                }
            }

            if(checkDrivers)
            {
                driverGeneAnnotator = new DriverGeneAnnotator(dbAccess, ensemblDataCache, config, cnDataLoader);
                driverGeneAnnotator.loadConfig(cmd);
                driverGeneAnnotator.setVisWriter(sampleAnalyser.getVisWriter());
            }
        }

        PerformanceCounter prefCounter = new PerformanceCounter("Total");

        int count = 0;
        for (final String sampleId : samplesList)
        {
            ++count;

            prefCounter.start();

            List<StructuralVariantData> svRecords = sampleDataFromFile ?
                    loadSampleSvDataFromFile(config.SvDataPath, sampleId, cmd) : dbAccess.readStructuralVariantData(sampleId);

            final List<SvVarData> svDataList = createSvData(svRecords);

            if(svDataList.isEmpty())
            {
                LOGGER.info("sample({}) has no passing SVs", sampleId);

                if(config.isSingleSample())
                {
                    sampleAnalyser.writeSampleWithNoSVs(sampleId);
                }

                continue;
            }

            if(config.hasMultipleSamples())
            {
                LOGGER.info("sample({}) processing {} SVs, completed({})", sampleId, svDataList.size(), count - 1);
            }

            cnDataLoader.loadSampleData(sampleId, svRecords);

            sampleAnalyser.setSampleSVs(sampleId, svDataList);

            if(ensemblDataCache != null)
            {
                ensemblDataCache.setSvGeneData(svDataList, applyPromotorDistance, selectiveGeneLoading);
            }

            sampleAnalyser.analyse();

            if(!sampleAnalyser.inValidState())
            {
                LOGGER.info("exiting after sample({}), in invalid state", sampleId);
                break;
            }

            if(checkDrivers || checkFusions)
            {
                // when matching RNA, allow all transcripts regardless of their viability for fusions
                fusionAnalyser.annotateTranscripts(svDataList, purgeInvalidTranscripts);
            }

            sampleAnalyser.annotate();

            if(checkDrivers)
            {
                driverGeneAnnotator.annotateSVs(sampleId, sampleAnalyser.getChrBreakendMap());
            }

            if(checkFusions)
            {
                fusionAnalyser.run(sampleId, svDataList, dbAccess, sampleAnalyser.getClusters(), sampleAnalyser.getChrBreakendMap());
            }

            sampleAnalyser.writeOutput(dbAccess);

            prefCounter.stop();

            if(config.MaxSamples > 0 && count >= config.MaxSamples)
            {
                LOGGER.info("exiting after max sample count {} reached", count);
                break;
            }
        }

        if(LOGGER.isDebugEnabled() || config.hasMultipleSamples())
        {
            prefCounter.logStats();
        }

        sampleAnalyser.close();

        if(fusionAnalyser != null)
            fusionAnalyser.close();

        if(driverGeneAnnotator != null)
            driverGeneAnnotator.close();

        if(config.isSingleSample())
        {
            try { version.write(config.OutputDataPath); } catch(IOException e) {}
        }

        LOGGER.info("SV analysis complete");
    }

    private static List<StructuralVariantData> loadSampleSvDataFromFile(final String samplePath, final String sampleId, final CommandLine cmd)
    {
        if(cmd.hasOption(VCF_FILE))
        {
            return loadSvDataFromVcf(cmd.getOptionValue(VCF_FILE));
        }
        else
        {
            return loadSvDataFromSvFile(sampleId, samplePath);
        }
    }

    private static List<SvVarData> createSvData(List<StructuralVariantData> svRecords)
    {
        List<SvVarData> svVarDataItems = Lists.newArrayList();

        for (final StructuralVariantData svRecord : svRecords)
        {
            if(svRecord.filter().isEmpty() || svRecord.filter().equals(PASS) || svRecord.filter().equals(INFERRED))
            {
                svVarDataItems.add(new SvVarData(svRecord));
            }
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
        options.addOption(VCF_FILE, true, "Path to the PURPLE structural variant VCF file");
        options.addOption(REF_GENOME_FILE, true, "Path to the indexed ref genome fasta file");

        // allow sub-components to add their specific config
        LinxConfig.addCmdLineArgs(options);
        FusionFinder.addCmdLineArgs(options);
        DriverGeneAnnotator.addCmdLineArgs(options);
        FusionDisruptionAnalyser.addCmdLineArgs(options);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
