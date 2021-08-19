package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.linx.LinxApplication.loadSampleSvDataFromFile;
import static com.hartwig.hmftools.linx.LinxConfig.CHECK_DRIVERS;
import static com.hartwig.hmftools.linx.LinxConfig.CHECK_FUSIONS;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.RG_VERSION;
import static com.hartwig.hmftools.linx.SvFileLoader.VCF_FILE;
import static com.hartwig.hmftools.linx.ext_compare.AmpliconCompare.AMPLICON_DATA_FILE;
import static com.hartwig.hmftools.linx.ext_compare.ChainFinderCompare.CHAIN_FINDER_DATA_DIR;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.analysis.CohortDataWriter;
import com.hartwig.hmftools.linx.analysis.SampleAnalyser;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.drivers.DriverGeneAnnotator;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class ExternalToolCompare
{
    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        if(!LinxConfig.validConfig(cmd))
        {
            LNX_LOGGER.error("exiting on invalid config");
            return;
        }

        LinxConfig config = new LinxConfig(cmd);

        final DatabaseAccess dbAccess = createDatabaseAccess(cmd);

        final List<String> samplesList = config.getSampleIds();

        boolean sampleDataFromFile = !config.PurpleDataPath.isEmpty() && cmd.hasOption(VCF_FILE);

        /*
        if(dbAccess == null)
        {
            LNX_LOGGER.error("no DB connection configured");
            return;
        }
        */

        if (samplesList.isEmpty())
        {
            LNX_LOGGER.error("samples must be specified");
            return;
        }

        LNX_LOGGER.info("running Linx external tool comparison for {} samples", samplesList.size());

        CohortDataWriter cohortDataWriter = new CohortDataWriter(config);

        SampleAnalyser sampleAnalyser = new SampleAnalyser(config, dbAccess, cohortDataWriter);

        CnDataLoader cnDataLoader = new CnDataLoader(config.PurpleDataPath, dbAccess);
        sampleAnalyser.setCnDataLoader(cnDataLoader);

        DriverGeneAnnotator driverGeneAnnotator = null;
        boolean checkDrivers = cmd.hasOption(CHECK_DRIVERS);

        FusionDisruptionAnalyser fusionAnalyser = null;

        ChainFinderCompare chainFinderCompare = cmd.hasOption(CHAIN_FINDER_DATA_DIR) ?
                new ChainFinderCompare(config.OutputDataPath, cmd) : null;

        AmpliconCompare ampliconCompare = cmd.hasOption(AMPLICON_DATA_FILE) ? new AmpliconCompare(config.OutputDataPath, cmd) : null;

        boolean selectiveGeneLoading = (samplesList.size() == 1 && !checkDrivers);

        final EnsemblDataCache ensemblDataCache = cmd.hasOption(GENE_TRANSCRIPTS_DIR) ?
                new EnsemblDataCache(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR), RG_VERSION) : null;

        if(ensemblDataCache != null)
        {
            ensemblDataCache.setRequiredData(true, false, false, true);

            if(!ensemblDataCache.load(selectiveGeneLoading))
            {
                LNX_LOGGER.error("Ensembl data cache load failed, exiting");
                return;
            }

            sampleAnalyser.setGeneCollection(ensemblDataCache);
            sampleAnalyser.getVisWriter().setGeneDataCache(ensemblDataCache);

            fusionAnalyser = new FusionDisruptionAnalyser(cmd, config, ensemblDataCache, cohortDataWriter);

            if(checkDrivers)
            {
                driverGeneAnnotator = new DriverGeneAnnotator(dbAccess, ensemblDataCache, config, cnDataLoader);
                driverGeneAnnotator.setVisWriter(sampleAnalyser.getVisWriter());
            }
        }

        int count = 0;
        for (final String sampleId : samplesList)
        {
            ++count;

            final List<StructuralVariantData> svRecords = sampleDataFromFile ?
                    loadSampleSvDataFromFile(config, sampleId, cmd) : dbAccess.readStructuralVariantData(sampleId);

            final List<SvVarData> svDataList = createSvData(svRecords);

            if(svDataList.isEmpty())
            {
                LNX_LOGGER.debug("sample({}) has no passing SVs", sampleId);
                continue;
            }

            if(config.hasMultipleSamples())
            {
                LNX_LOGGER.info("sample({}) processing {} SVs, completed({})", sampleId, svDataList.size(), count - 1);
            }

            sampleAnalyser.setSampleId(sampleId);
            cnDataLoader.loadSampleData(sampleId, svRecords);
            sampleAnalyser.setSampleSVs(svDataList);

            if(ensemblDataCache != null)
            {
                sampleAnalyser.setSvGeneData(svDataList, ensemblDataCache, false, selectiveGeneLoading);
            }

            sampleAnalyser.analyse();

            if(!sampleAnalyser.inValidState())
            {
                LNX_LOGGER.info("exiting after sample({}), in invalid state", sampleId);
                break;
            }

            if(checkDrivers)
            {
                fusionAnalyser.annotateTranscripts(svDataList, false);
            }

            sampleAnalyser.annotate();

            if(checkDrivers)
            {
                driverGeneAnnotator.annotateSVs(sampleId, sampleAnalyser.getChrBreakendMap());
            }

            sampleAnalyser.writeOutput(dbAccess);

            if(chainFinderCompare != null)
            {
                chainFinderCompare.processSample(sampleId, svDataList, sampleAnalyser.getClusters(), sampleAnalyser.getChrBreakendMap());
            }

            if(ampliconCompare != null)
            {
                ampliconCompare.processSample(sampleId, sampleAnalyser.getChrBreakendMap());
            }
        }

        sampleAnalyser.close();

        if(driverGeneAnnotator != null)
            driverGeneAnnotator.close();

        if(chainFinderCompare != null)
            chainFinderCompare.close();

        if(ampliconCompare != null)
            ampliconCompare.close();

        LNX_LOGGER.info("external tool comparison complete");
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

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        addDatabaseCmdLineArgs(options);
        options.addOption(CHECK_DRIVERS, false, "Check SVs against drivers catalog");
        options.addOption(CHECK_FUSIONS, false, "Run fusion detection");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Optional: Ensembl data cache directory");
        options.addOption(VCF_FILE, true, "Path to the PURPLE structural variant VCF file");

        // allow sub-components to add their specific config
        LinxConfig.addCmdLineArgs(options);
        ChainFinderCompare.addCmdLineArgs(options);
        AmpliconCompare.addCmdLineArgs(options);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
