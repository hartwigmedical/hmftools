package com.hartwig.hmftools.linx;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.REF_GENOME_VERSION;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.linx.fusion.FusionConfig;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.fusion.FusionResources;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LinxApplication
{
    private static final String FILTER_QC_PASS = "filter_qc_pass";

    public LinxApplication(final CommandLine cmd)
    {
        final VersionInfo version = new VersionInfo("linx.version");
        LNX_LOGGER.info("LINX version: {}", version.version());

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

        long startTime = System.currentTimeMillis();

        final DatabaseAccess dbAccess = createDatabaseAccess(cmd);

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
            System.exit(1);
        }

        LNX_LOGGER.info("running {} analysis for {}",
                config.IsGermline ? "germline SV" : "SV",
                config.hasMultipleSamples() ? String.format("%d samples", samplesList.size()) : samplesList.get(0));

        FusionResources fusionResources = new FusionResources(cmd);

        if(config.RunFusions && !fusionResources.knownFusionCache().hasValidData())
        {
            LNX_LOGGER.info("invalid fusion config, exiting");
            System.exit(1);
        }

        if(!cmd.hasOption(ENSEMBL_DATA_DIR))
        {
            if(config.RunFusions || config.RunDrivers || config.IsGermline)
            {
                LNX_LOGGER.info("missing Ensembl data cache, exiting");
                System.exit(1);
            }
        }

        final EnsemblDataCache ensemblDataCache = cmd.hasOption(ENSEMBL_DATA_DIR) ?
                new EnsemblDataCache(cmd.getOptionValue(ENSEMBL_DATA_DIR), REF_GENOME_VERSION) : null;

        if(ensemblDataCache != null)
        {
            boolean reqProteinDomains = config.RunFusions;
            boolean reqSplicePositions = config.RunFusions;
            boolean canonicalOnly = config.IsGermline;

            ensemblDataCache.setRequiredData(true, reqProteinDomains, reqSplicePositions, canonicalOnly);

            boolean ensemblLoadOk = false;

            if(!config.RestrictedGeneIds.isEmpty())
            {
                ensemblDataCache.setRestrictedGeneIdList(config.RestrictedGeneIds);
                ensemblLoadOk = ensemblDataCache.load(false);
            }
            else if(!config.RunFusions && config.RunDrivers && !config.HomDisAllGenes)
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

            if(!ensemblLoadOk)
            {
                LNX_LOGGER.error("Ensembl data cache load failed, exiting");
                System.exit(1);
            }

            if(config.RunDrivers)
                ensemblDataCache.createGeneNameIdMap();
        }

        CohortDataWriter cohortDataWriter = new CohortDataWriter(config, ensemblDataCache);

        SvAnnotators svAnnotators = new SvAnnotators(config, ensemblDataCache, dbAccess);

        List<SampleAnalyser> sampleAnalysers = Lists.newArrayList();

        if(config.Threads > 1)
        {
            List<List<String>> saSampleLists = Lists.newArrayList();

            int threads = min(config.Threads, samplesList.size());

            for(int i = 0; i < threads; ++i)
            {
                sampleAnalysers.add(new SampleAnalyser(
                        i, config, dbAccess, svAnnotators, ensemblDataCache, fusionResources, cohortDataWriter));

                saSampleLists.add(Lists.newArrayList());
            }

            int saIndex = 0;
            for (final String sampleId : samplesList)
            {
                saSampleLists.get(saIndex).add(sampleId);
                ++saIndex;

                if(saIndex >= threads)
                    saIndex = 0;
            }

            for(int i = 0; i < threads; ++i)
            {
                sampleAnalysers.get(i).setSampleIds(saSampleLists.get(i));
            }

            final List<Callable> callableList = sampleAnalysers.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, callableList.size());
        }
        else
        {
            SampleAnalyser sampleAnalyser = new SampleAnalyser(
                    0, config, dbAccess, svAnnotators, ensemblDataCache, fusionResources, cohortDataWriter);

            sampleAnalysers.add(sampleAnalyser);
            sampleAnalyser.setSampleIds(config.getSampleIds());
            sampleAnalyser.processSamples();
        }

        svAnnotators.close();
        cohortDataWriter.close();

        if(config.hasMultipleSamples())
        {
            // combine and log performance counters
            Map<String,PerformanceCounter> combinedPerfCounters = sampleAnalysers.get(0).getPerfCounters();

            for(int i = 1; i < sampleAnalysers.size(); ++i)
            {
                Map<String,PerformanceCounter> saPerfCounters = sampleAnalysers.get(i).getPerfCounters();

                for(Map.Entry<String,PerformanceCounter> entry : combinedPerfCounters.entrySet())
                {
                    PerformanceCounter combinedPc = entry.getValue();
                    PerformanceCounter saPc = saPerfCounters.get(entry.getKey());

                    if(combinedPc != null)
                        combinedPc.merge(saPc);
                }
            }

            for(Map.Entry<String,PerformanceCounter> entry : combinedPerfCounters.entrySet())
            {
                entry.getValue().logStats();
            }
        }

        if(config.isSingleSample())
        {
            try { version.write(config.OutputDataPath); } catch(IOException e) {}
        }

        if(config.isSingleSample())
        {
            LNX_LOGGER.info("Linx complete for {}", samplesList.get(0));
        }
        else
        {
            double runTime = (System.currentTimeMillis() - startTime) / 1000.0;

            LNX_LOGGER.info("SV analysis complete for {} samples, run time({})s",
                    samplesList.size(), String.format("%.3f", runTime));
        }
    }

    private static List<String> getStructuralVariantSamplesList(@NotNull DatabaseAccess dbAccess, boolean filterQCPassOnly)
    {
        final List<String> sampleIds = filterQCPassOnly ? dbAccess.readPurpleSampleListPassingQC(MIN_SAMPLE_PURITY) : dbAccess.readPurpleSampleList();

        if(!sampleIds.isEmpty())
            return sampleIds;

        return dbAccess.readStructuralVariantSampleList("");
    }

    public static void main(@NotNull final String[] args)
    {
        final Options options = createBasicOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            new LinxApplication(cmd);
        }
        catch(ParseException e)
        {
            LNX_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("LinxApplication", options);
            System.exit(1);
        }

    }

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        addDatabaseCmdLineArgs(options);
        addEnsemblDir(options);
        options.addOption(FILTER_QC_PASS, false, "Optional: If present will filter out QC-fail sample");

        // allow sub-components to add their specific config
        LinxConfig.addCmdLineArgs(options);
        addKnownFusionFileOption(options);
        FusionConfig.addCmdLineArgs(options);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
