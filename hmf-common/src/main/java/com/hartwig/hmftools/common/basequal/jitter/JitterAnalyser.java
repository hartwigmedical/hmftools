package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.basequal.jitter.RefGenomeMicrosatellite.filterMicrosatellites;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.bam.BamSlicerFilter;
import com.hartwig.hmftools.common.utils.r.RExecutor;

import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class JitterAnalyser
{
    private final JitterAnalyserConfig mConfig;
    private final Logger mLogger;
    private final BamSlicerFilter mBamSlicerFilter;
    private final SampleReadProcessor mSampleReadProcessor;

    public JitterAnalyser(final JitterAnalyserConfig config, final Logger logger)
    {
        mConfig = config;
        mLogger = logger;

        mBamSlicerFilter = new BamSlicerFilter(config.MinMappingQuality, false, false, false);

        List<RefGenomeMicrosatellite> refGenomeMicrosatellites = loadRefGenomeMicrosatellites();
        mSampleReadProcessor = new SampleReadProcessor(refGenomeMicrosatellites);
    }

    public void processBam() throws InterruptedException
    {
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        SampleBamProcessor sampleBamProcessor = new SampleBamProcessor(mConfig, mBamSlicerFilter, mSampleReadProcessor);
        sampleBamProcessor.queryBam(mConfig, executorService);
    }

    public void processRead(final SAMRecord read)
    {
        mSampleReadProcessor.processRead(read);
    }

    public void writeAnalysisOutput() throws IOException, InterruptedException
    {
        Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers = mSampleReadProcessor.getMicrosatelliteSiteAnalysers();

        // now write out all the repeat stats
        MicrosatelliteSiteFile.write(MicrosatelliteSiteFile.generateFilename(mConfig.OutputDir, mConfig.SampleId), microsatelliteSiteAnalysers);

        final String statsTableFile = JitterCountsTableFile.generateFilename(mConfig.OutputDir, mConfig.SampleId);
        writeMicrosatelliteStatsTable(microsatelliteSiteAnalysers, statsTableFile);

        // draw a chart of the 9 ms profiles
        if(mConfig.WritePlots)
            drawMicrosatelliteCharts(mConfig.OutputDir, mConfig.SampleId, statsTableFile);

        // now perform the fitting
        List<JitterModelParams> jitterModelParamsList = fitJitterModels(microsatelliteSiteAnalysers);

        JitterModelParamsFile.write(JitterModelParamsFile.generateFilename(mConfig.OutputDir, mConfig.SampleId), jitterModelParamsList);
    }

    private List<RefGenomeMicrosatellite> loadRefGenomeMicrosatellites()
    {
        List<RefGenomeMicrosatellite> refGenomeMicrosatellites = RefGenomeMicrosatelliteFile.read(mConfig.RefGenomeMsiFile);

        mLogger.info("loaded {} microsatellites regions", refGenomeMicrosatellites.size());

        filterSpecificRegions(refGenomeMicrosatellites);
        refGenomeMicrosatellites = filterMicrosatellites(refGenomeMicrosatellites, mConfig.MaxSitesPerType);

        return refGenomeMicrosatellites;
    }

    private static void writeMicrosatelliteStatsTable(@NotNull final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers,
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

    private static void filterSpecificRegions(List<RefGenomeMicrosatellite> refGenomeMicrosatellites)
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
        int result = RExecutor.executeFromClasspath("basequal/msi_jitter_plot.R", outputDir, sampleId, statsTableFile);
        if(result != 0)
        {
            throw new IOException("R execution failed. Unable to complete segmentation.");
        }
    }

    public BamSlicerFilter bamSlicerFilter()
    {
        return mBamSlicerFilter;
    }
}
