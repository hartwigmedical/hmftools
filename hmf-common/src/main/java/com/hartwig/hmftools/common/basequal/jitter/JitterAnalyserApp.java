package com.hartwig.hmftools.common.basequal.jitter;

import static java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME;

import static com.hartwig.hmftools.common.basequal.jitter.RefGenomeMicrosatellite.filterMicrosatellites;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;

import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class JitterAnalyserApp
{
    public static final Logger sLogger = LogManager.getLogger(JitterAnalyserApp.class);

    private final JitterAnalyserConfig mConfig;

    public JitterAnalyserApp(final ConfigBuilder configBuilder) throws ParseException
    {
        mConfig = new JitterAnalyserConfig(configBuilder);
    }

    public int run(final String... args) throws InterruptedException, IOException
    {
        Instant start = Instant.now();

        VersionInfo versionInfo = new VersionInfo("errorprofile.version");

        sLogger.info("ErrorProfile version: {}", versionInfo.version());

        sLogger.debug("build timestamp: {}, run args: {}",
                versionInfo.buildTime().format(ISO_ZONED_DATE_TIME), String.join(" ", args));

        if(!mConfig.isValid())
        {
            sLogger.error(" invalid config, exiting");
            return 1;
        }

        List<RefGenomeMicrosatellite> refGenomeMicrosatellites = RefGenomeMicrosatelliteFile.read(mConfig.RefGenomeMicrosatelliteFile);

        sLogger.info("loaded {} microsatellites regions", refGenomeMicrosatellites.size());

        filterSpecificRegions(refGenomeMicrosatellites);
        refGenomeMicrosatellites = filterMicrosatellites(refGenomeMicrosatellites, mConfig.MaxSitesPerType);

        SampleBamProcessor sampleBamProcessor = new SampleBamProcessor(refGenomeMicrosatellites);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        sampleBamProcessor.queryBam(mConfig, executorService);

        // now write out all the repeat stats
        MicrosatelliteSiteFile.write(MicrosatelliteSiteFile.generateFilename(mConfig.OutputDir, mConfig.SampleId),
                sampleBamProcessor.getMicrosatelliteSiteAnalysers());

        final String statsTableFile = JitterCountsTableFile.generateFilename(mConfig.OutputDir, mConfig.SampleId);
        writeMicrosatelliteStatsTable(sampleBamProcessor.getMicrosatelliteSiteAnalysers(), statsTableFile);

        // draw a chart of the 9 ms profiles
        drawMicrosatelliteCharts(mConfig.OutputDir, mConfig.SampleId, statsTableFile);

        // now perform the fitting
        List<JitterModelParams> jitterModelParamsList = fitJitterModels(sampleBamProcessor.getMicrosatelliteSiteAnalysers());

        JitterModelParamsFile.write(JitterModelParamsFile.generateFilename(mConfig.OutputDir, mConfig.SampleId), jitterModelParamsList);

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        sLogger.info("run complete. Time taken: {}m {}s", seconds / 60, seconds % 60);

        return 0;
    }

    private void writeMicrosatelliteStatsTable(@NotNull final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers,
            final String filename)
    {
        // write two tables, one with real variant filter, one without

        List<JitterCountsTable> msStatsTables = new ArrayList<>();

        for(MicrosatelliteSelector s : createMicrosatelliteSelectorsForCharts())
        {
            JitterCountsTable msStatsTable = JitterCountsTable.summariseFrom(s.unitName(),
                microsatelliteSiteAnalysers.stream().filter(s::select).collect(Collectors.toList()));
            msStatsTables.add(msStatsTable);
        }

        JitterCountsTableFile.write(filename, msStatsTables);
    }

    private static List<MicrosatelliteSelector> createMicrosatelliteSelectorsForCharts()
    {
        // create nine summary / pivot table
        // {A/T, C/G, AT/TA, AG/GA/CT/TC, AC/CA/GT/TG, CG/GC, any 3 base, any 4 base, any 5 base}
        List<MicrosatelliteSelector> selectors = new ArrayList<>();
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("A", "T")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("C", "G")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("AT", "TA")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("AG", "GA", "CT", "TC")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("AC", "CA", "GT", "TG")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("CG", "GC")));
        selectors.add(MicrosatelliteSelector.fromUnitLength(3));
        selectors.add(MicrosatelliteSelector.fromUnitLength(4));
        selectors.add(MicrosatelliteSelector.fromUnitLength(5));
        //selectors.add(MicrosatelliteSelector.fromUnitLengthRange(3, 5));

        return selectors;
    }

    private static List<JitterModelParams> fitJitterModels(
            @NotNull final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers)
    {
        // create nine summary / pivot table
        // {A/T, C/G, AT/TA, AG/GA/CT/TC, AC/CA/GT/TG, CG/GC, any 3 base, any 4 base, any 5 base}
        List<MicrosatelliteSelector> selectors = new ArrayList<>();
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("A", "T")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("C", "G")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("AT", "TA")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("AG", "GA", "CT", "TC")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("AC", "CA", "GT", "TG")));
        selectors.add(MicrosatelliteSelector.fromUnits(List.of("CG", "GC")));
        selectors.add(MicrosatelliteSelector.fromUnitLengthRange(3, 5));

        List<JitterModelParams> fittedParams = new ArrayList<>();

        for(MicrosatelliteSelector selector : selectors)
        {
            JitterCountsTable msStatsTable = JitterCountsTable.summariseFrom(selector.unitName(),
                microsatelliteSiteAnalysers.stream().filter(selector::select).collect(Collectors.toList()));
            JitterModelFitter fitter = new JitterModelFitter(msStatsTable);
            fitter.performFit();
            fittedParams.add(fitter.getJitterModelParams());
        }

        return fittedParams;
    }

    private void filterSpecificRegions(List<RefGenomeMicrosatellite> refGenomeMicrosatellites)
    {
        /*
        if(!mConfig.SpecificRegions.isEmpty())
        {
            sLogger.info(mConfig.SpecificRegions);

            refGenomeMicrosatellites.removeIf(refGenomeHomopolymer -> mConfig.SpecificRegions.stream()
                    .noneMatch(o -> o.overlaps(refGenomeHomopolymer.genomeRegion)));
        }
         */

        // we want to get around 1000 sites for each repeat context
    }

    private static void drawMicrosatelliteCharts(final String outputDir, final String sampleId, final String statsTableFile)
            throws IOException, InterruptedException
    {
        int result = RExecutor.executeFromClasspath("basequal/microsatellite/ms_plot.R", outputDir, sampleId, statsTableFile);
        if(result != 0)
        {
            throw new IOException("R execution failed. Unable to complete segmentation.");
        }
    }

    public static void main(final String... args) throws InterruptedException, IOException, ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("ErrorProfile");
        JitterAnalyserConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        // set all thread exception handler
        // if we do not do this, exception thrown in other threads will not be handled and results
        // in the program hanging
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            sLogger.error("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        System.exit(new JitterAnalyserApp(configBuilder).run(args));
    }
}
