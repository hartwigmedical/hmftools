package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.linx.LinxConfig.CHECK_DRIVERS;
import static com.hartwig.hmftools.linx.LinxConfig.CHECK_FUSIONS;
import static com.hartwig.hmftools.linx.LinxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.linx.LinxConfig.LOG_VERBOSE;
import static com.hartwig.hmftools.linx.LinxConfig.RG_VERSION;
import static com.hartwig.hmftools.linx.SvFileLoader.VCF_FILE;
import static com.hartwig.hmftools.linx.SvFileLoader.loadSvDataFromGermlineVcf;
import static com.hartwig.hmftools.linx.SvFileLoader.loadSvDataFromVcf;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.linx.analysis.SampleAnalyser;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.drivers.DriverGeneAnnotator;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.fusion.FusionFinder;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class LinxApplication
{
    private static final String FILTER_QC_PASS = "filter_qc_pass";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final VersionInfo version = new VersionInfo("linx.version");
        LNX_LOGGER.info("LINX version: {}", version.version());

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

        if(!LinxConfig.validConfig(cmd))
        {
            LNX_LOGGER.error("exiting on invalid config");
            return;
        }

        LinxConfig config = new LinxConfig(cmd);

        if(!checkCreateOutputDir(config.OutputDataPath))
        {
            LNX_LOGGER.error("failed to create output directory({})", config.OutputDataPath);
            return;
        }

        final DatabaseAccess dbAccess = createDatabaseAccess(cmd);

        boolean sampleDataFromFile = (!config.PurpleDataPath.isEmpty() && config.SvVcfFile != null) || config.IsGermline;

        List<String> samplesList = config.getSampleIds();

        if(dbAccess == null && !config.hasValidSampleDataSource(cmd))
            return;

        if (samplesList.isEmpty())
        {
            boolean filterQCPassOnly = cmd.hasOption(FILTER_QC_PASS);
            samplesList = getStructuralVariantSamplesList(dbAccess, filterQCPassOnly);

            LNX_LOGGER.info("retrieved {} samples {}", samplesList.size(), filterQCPassOnly ? "QC-pass-filtered" : "");

            config.setSampleIds(samplesList);
        }

        if (samplesList.isEmpty())
        {
            LNX_LOGGER.info("not samples loaded, exiting");
            return;
        }

        LNX_LOGGER.info("running SV analysis for {}",
                config.hasMultipleSamples() ? String.format("%d samples", samplesList.size()) : samplesList.get(0));

        SampleAnalyser sampleAnalyser = new SampleAnalyser(config, dbAccess);

        CnDataLoader cnDataLoader = new CnDataLoader(config.PurpleDataPath, dbAccess);

        if(config.hasMultipleSamples())
            cnDataLoader.establishCaches();

        sampleAnalyser.setCnDataLoader(cnDataLoader);

        DriverGeneAnnotator driverGeneAnnotator = null;
        boolean checkDrivers = cmd.hasOption(CHECK_DRIVERS) && config.DriverGenes != null;

        FusionDisruptionAnalyser fusionAnalyser = null;
        boolean checkFusions = cmd.hasOption(CHECK_FUSIONS);

        boolean breakendGeneLoading = (samplesList.size() == 1 && !checkDrivers) && config.RestrictedGeneIds.isEmpty();
        boolean applyPromotorDistance = checkFusions;
        boolean purgeInvalidTranscripts = true;

        final EnsemblDataCache ensemblDataCache = cmd.hasOption(GENE_TRANSCRIPTS_DIR) ?
                new EnsemblDataCache(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR), RG_VERSION) : null;

        if(ensemblDataCache != null)
        {
            ensemblDataCache.setRequiredData(true, checkFusions, checkFusions, !checkFusions);

            boolean ensemblLoadOk = false;

            if(!breakendGeneLoading)
            {
                if(!config.RestrictedGeneIds.isEmpty())
                {
                    ensemblDataCache.setRestrictedGeneIdList(config.RestrictedGeneIds);
                    ensemblLoadOk = ensemblDataCache.load(false);
                }
                else if(!checkFusions && checkDrivers)
                {
                    // only load transcripts for the driver gene panel
                    ensemblLoadOk = ensemblDataCache.load(true);

                    if(ensemblLoadOk)
                    {
                        final List<String> geneIds = config.DriverGenes.stream()
                                .map(x -> ensemblDataCache.getGeneDataByName(x.gene()))
                                .filter(x -> x != null)
                                .map(x -> x.GeneId).collect(Collectors.toList());

                        ensemblLoadOk &= ensemblDataCache.loadTranscriptData(geneIds);
                    }
                }
                else
                {
                    ensemblLoadOk = ensemblDataCache.load(false);
                }
            }
            else
            {
                ensemblLoadOk = ensemblDataCache.load(true);
            }

            if(!ensemblLoadOk)
            {
                LNX_LOGGER.error("Ensembl data cache load failed, exiting");
                return;
            }

            sampleAnalyser.setGeneCollection(ensemblDataCache);
            sampleAnalyser.getVisWriter().setGeneDataCache(ensemblDataCache);

            // always initialise since is used for transcript evaluation
            fusionAnalyser = new FusionDisruptionAnalyser(cmd, config, ensemblDataCache, sampleAnalyser.getVisWriter());

            if(!fusionAnalyser.validState())
                return;

            if(checkFusions)
            {
                purgeInvalidTranscripts = !fusionAnalyser.hasRnaSampleData();

                if(fusionAnalyser.hasRnaSampleData() && samplesList.size() > 1)
                {
                    samplesList.clear();
                    samplesList.addAll(fusionAnalyser.getRnaSampleIds());

                    LNX_LOGGER.info("running {} sample based on RNA fusion input", samplesList.size());
                }
            }

            if(checkDrivers)
            {
                driverGeneAnnotator = new DriverGeneAnnotator(dbAccess, ensemblDataCache, config, cnDataLoader);
                driverGeneAnnotator.setVisWriter(sampleAnalyser.getVisWriter());
            }
        }

        PerformanceCounter prefCounter = new PerformanceCounter("Total");

        int count = 0;
        for (final String sampleId : samplesList)
        {
            ++count;

            prefCounter.start();

            final List<StructuralVariantData> svRecords = sampleDataFromFile ?
                    loadSampleSvDataFromFile(config, sampleId, cmd) : dbAccess.readStructuralVariantData(sampleId);

            final List<SvVarData> svDataList = createSvData(svRecords, config);

            sampleAnalyser.setSampleId(sampleId);

            if(svDataList.isEmpty())
            {
                LNX_LOGGER.info("sample({}) has no passing SVs", sampleId);

                if(config.isSingleSample())
                    sampleAnalyser.writeSampleWithNoSVs();

                continue;
            }

            if(config.hasMultipleSamples())
            {
                LNX_LOGGER.info("sample({}) processing {} SVs, completed({})", sampleId, svDataList.size(), count - 1);
            }

            if(!config.IsGermline)
                cnDataLoader.loadSampleData(sampleId, svRecords);

            sampleAnalyser.setSampleSVs(svDataList);

            if(ensemblDataCache != null)
            {
                sampleAnalyser.setSvGeneData(svDataList, ensemblDataCache, applyPromotorDistance, breakendGeneLoading);
            }

            sampleAnalyser.analyse();

            if(!sampleAnalyser.inValidState())
            {
                LNX_LOGGER.info("exiting after sample({}), in invalid state", sampleId);
                break;
            }

            if(fusionAnalyser != null)
            {
                // when matching RNA, allow all transcripts regardless of their viability for fusions
                fusionAnalyser.annotateTranscripts(svDataList, purgeInvalidTranscripts);
            }

            sampleAnalyser.annotate();

            if(checkDrivers)
            {
                driverGeneAnnotator.annotateSVs(sampleId, sampleAnalyser.getChrBreakendMap());
            }

            if(checkFusions || config.IsGermline)
            {
                fusionAnalyser.run(sampleId, svDataList, dbAccess, sampleAnalyser.getClusters(), sampleAnalyser.getChrBreakendMap());
            }

            sampleAnalyser.writeOutput(dbAccess);

            prefCounter.stop();
        }

        if(LNX_LOGGER.isDebugEnabled() || config.hasMultipleSamples())
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

        LNX_LOGGER.info("SV analysis complete for {}",
                config.hasMultipleSamples() ? String.format("%d samples", samplesList.size()) : samplesList.get(0));
    }

    public static List<StructuralVariantData> loadSampleSvDataFromFile(
            final LinxConfig config, final String sampleId, final CommandLine cmd)
    {
        String vcfFile = config.SvVcfFile.contains("*") ? config.SvVcfFile.replaceAll("\\*", sampleId) : config.SvVcfFile;

        if(config.IsGermline)
            return loadSvDataFromGermlineVcf(vcfFile);
        else
            return loadSvDataFromVcf(vcfFile);

        // no longer supported
        // return loadSvDataFromSvFile(sampleId, config.SampleDataPath);
    }

    private static List<SvVarData> createSvData(final List<StructuralVariantData> svRecords, final LinxConfig config)
    {
        List<SvVarData> svVarDataItems = Lists.newArrayList();

        for (final StructuralVariantData svRecord : svRecords)
        {
            final String filter = svRecord.filter();

            if(filter.isEmpty() || filter.equals(PASS) || filter.equals(INFERRED) || (config.IsGermline && filter.equals(PON_FILTER_PON)))
            {
                svVarDataItems.add(new SvVarData(svRecord));
            }
        }

        return svVarDataItems;
    }


    private static List<String> getStructuralVariantSamplesList(@NotNull DatabaseAccess dbAccess, boolean filterQCPassOnly)
    {
        final List<String> sampleIds = filterQCPassOnly ? dbAccess.readPurpleSampleListPassingQC(MIN_SAMPLE_PURITY) : dbAccess.readPurpleSampleList();

        if(!sampleIds.isEmpty())
            return sampleIds;

        return dbAccess.readStructuralVariantSampleList("");
    }

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        addDatabaseCmdLineArgs(options);
        options.addOption(CHECK_DRIVERS, false, "Check SVs against drivers catalog");
        options.addOption(CHECK_FUSIONS, false, "Run fusion detection");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Optional: Ensembl data cache directory");
        options.addOption(FILTER_QC_PASS, false, "Optional: If present will filter out QC-fail sample");
        options.addOption(VCF_FILE, true, "Path to the PURPLE structural variant VCF file");

        // allow sub-components to add their specific config
        LinxConfig.addCmdLineArgs(options);
        FusionFinder.addCmdLineArgs(options);
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
