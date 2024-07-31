package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.DUAL_BASE_1;
import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.DUAL_BASE_2;
import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.DUAL_BASE_3;
import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.DUAL_BASE_4;
import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.SINGLE_BASE_1;
import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.SINGLE_BASE_2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bam.BamSlicerFilter;
import com.hartwig.hmftools.common.utils.r.RExecutor;

import org.apache.logging.log4j.Logger;

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

    public BamSlicerFilter bamSlicerFilter()
    {
        return mBamSlicerFilter;
    }
    public SampleReadProcessor sampleReadProcessor() { return mSampleReadProcessor; }

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
        writeMicrosatelliteStatsTable(microsatelliteSiteAnalysers, statsTableFile, mConfig);

        // draw a chart of the 9 ms profiles
        if(mConfig.WritePlots)
            drawMicrosatelliteCharts(mConfig.OutputDir, mConfig.SampleId, statsTableFile);

        // now perform the fitting
        List<JitterModelParams> jitterModelParamsList = fitJitterModels(microsatelliteSiteAnalysers, mConfig.MaxSingleSiteAltContribution);

        JitterModelParamsFile.write(JitterModelParamsFile.generateFilename(mConfig.OutputDir, mConfig.SampleId), jitterModelParamsList);
    }

    private List<RefGenomeMicrosatellite> loadRefGenomeMicrosatellites()
    {
        List<RefGenomeMicrosatellite> refGenomeMicrosatellites = RefGenomeMicrosatelliteFile.read(mConfig.RefGenomeMsiFile);

        mLogger.info("loaded {} microsatellites regions", refGenomeMicrosatellites.size());

        filterSpecificRegions(refGenomeMicrosatellites);
        // refGenomeMicrosatellites = filterMicrosatellites(refGenomeMicrosatellites, mConfig.MaxSitesPerType);

        return refGenomeMicrosatellites;
    }

    private static void writeMicrosatelliteStatsTable(
            final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers, final String filename,
            final JitterAnalyserConfig config)
    {
        // write two tables, one with real variant filter, one without

        List<JitterCountsTable> msStatsTables = new ArrayList<>();

        for(MicrosatelliteSelector s : createMicrosatelliteSelectorsForCharts())
        {
            JitterCountsTable msStatsTable = JitterCountsTable.summariseFrom(
                    s.unitName(), config.MaxSingleSiteAltContribution,
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
        selectors.add(MicrosatelliteSelector.fromUnits(SINGLE_BASE_1));
        selectors.add(MicrosatelliteSelector.fromUnits(SINGLE_BASE_2));
        selectors.add(MicrosatelliteSelector.fromUnits(DUAL_BASE_1));
        selectors.add(MicrosatelliteSelector.fromUnits(DUAL_BASE_2));
        selectors.add(MicrosatelliteSelector.fromUnits(DUAL_BASE_3));
        selectors.add(MicrosatelliteSelector.fromUnits(DUAL_BASE_4));
        selectors.add(MicrosatelliteSelector.fromUnitLength(3));
        selectors.add(MicrosatelliteSelector.fromUnitLength(4));
        selectors.add(MicrosatelliteSelector.fromUnitLength(5));
        //selectors.add(MicrosatelliteSelector.fromUnitLengthRange(3, 5));

        return selectors;
    }

    private static List<JitterModelParams> fitJitterModels(
            final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers, final double maxSingleSiteAltContribution)
    {
        // create nine summary / pivot table
        // {A/T, C/G, AT/TA, AG/GA/CT/TC, AC/CA/GT/TG, CG/GC, any 3 base, any 4 base, any 5 base}
        List<MicrosatelliteSelector> selectors = new ArrayList<>();
        selectors.add(MicrosatelliteSelector.fromUnits(SINGLE_BASE_1));
        selectors.add(MicrosatelliteSelector.fromUnits(SINGLE_BASE_2));
        selectors.add(MicrosatelliteSelector.fromUnits(DUAL_BASE_1));
        selectors.add(MicrosatelliteSelector.fromUnits(DUAL_BASE_2));
        selectors.add(MicrosatelliteSelector.fromUnits(DUAL_BASE_3));
        selectors.add(MicrosatelliteSelector.fromUnits(DUAL_BASE_4));
        selectors.add(MicrosatelliteSelector.fromUnitLengthRange(3, 5));

        List<JitterModelParams> fittedParams = new ArrayList<>();

        for(MicrosatelliteSelector selector : selectors)
        {
            JitterCountsTable msStatsTable = JitterCountsTable.summariseFrom(
                    selector.unitName(), maxSingleSiteAltContribution,
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
}
