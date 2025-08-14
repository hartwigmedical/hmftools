package com.hartwig.hmftools.redux.jitter;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.jitter.JitterAnalyserConstants.DUAL_BASE_1;
import static com.hartwig.hmftools.redux.jitter.JitterAnalyserConstants.DUAL_BASE_2;
import static com.hartwig.hmftools.redux.jitter.JitterAnalyserConstants.DUAL_BASE_3;
import static com.hartwig.hmftools.redux.jitter.JitterAnalyserConstants.DUAL_BASE_4;
import static com.hartwig.hmftools.redux.jitter.JitterAnalyserConstants.SINGLE_BASE_1;
import static com.hartwig.hmftools.redux.jitter.JitterAnalyserConstants.SINGLE_BASE_2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamSlicerFilter;
import com.hartwig.hmftools.common.redux.JitterCountsTable;
import com.hartwig.hmftools.common.redux.JitterCountsTableFile;
import com.hartwig.hmftools.common.redux.JitterModelParams;
import com.hartwig.hmftools.common.redux.JitterModelParamsFile;
import com.hartwig.hmftools.common.redux.JitterTableRow;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.utils.RExecutor;

import htsjdk.samtools.SAMRecord;

public class JitterAnalyser
{
    private final JitterAnalyserConfig mConfig;
    private final BamSlicerFilter mBamSlicerFilter;
    private final SampleReadProcessor mSampleReadProcessor;

    private EnumSet<ConsensusType> mConsensusTypes;

    public JitterAnalyser(final JitterAnalyserConfig config)
    {
        mConfig = config;

        mBamSlicerFilter = new BamSlicerFilter(config.MinMappingQuality, false, false, false);

        List<RefGenomeMicrosatellite> refGenomeMicrosatellites = loadRefGenomeMicrosatellites();
        ConsensusMarker consensusMarker = ConsensusMarker.create(config);
        mSampleReadProcessor = new SampleReadProcessor(config, refGenomeMicrosatellites, consensusMarker);

        mConsensusTypes = null;
    }

    public BamSlicerFilter bamSlicerFilter()
    {
        return mBamSlicerFilter;
    }

    public void processRead(final SAMRecord read)
    {
        mSampleReadProcessor.processRead(read);
    }

    private EnumSet<ConsensusType> consensusTypes()
    {
        if(mConsensusTypes != null)
            return mConsensusTypes;

        mConsensusTypes = JitterAnalyserConfig.consensusTypes(mConfig);
        return mConsensusTypes;
    }

    public void writeAnalysisOutput() throws IOException, InterruptedException
    {
        Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers = mSampleReadProcessor.getMicrosatelliteSiteAnalysers();

        // now write out all the repeat stats
        if(mConfig.WriteSiteFile)
            MicrosatelliteSiteFile.write(MicrosatelliteSiteFile.generateFilename(mConfig.OutputDir, mConfig.SampleId), microsatelliteSiteAnalysers, consensusTypes());

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
        RD_LOGGER.info("loaded {} microsatellites regions", refGenomeMicrosatellites.size());
        return refGenomeMicrosatellites;
    }

    private void writeMicrosatelliteStatsTable(
            final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers, final String filename)
    {
        // write two tables, one with real variant filter, one without

        List<JitterCountsTable> msStatsTables = new ArrayList<>();

        for(ConsensusType consensusType : consensusTypes())
        {
            for(MicrosatelliteSelector s : createMicrosatelliteSelectorsForCharts())
            {
                JitterCountsTable msStatsTable = buildJitterCountsTable(
                        s.unitName(), consensusType, microsatelliteSiteAnalysers.stream().filter(s::select).collect(Collectors.toList()));
                msStatsTables.add(msStatsTable);
            }
        }

        JitterCountsTableFile.write(filename, msStatsTables);
    }

    static JitterCountsTable buildJitterCountsTable(
            final String repeatUnit, final ConsensusType consensusType,
            final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers)
    {
        // In order to filter out outliers, we perform the stats summation in a loop
        // We create a table of all read stats, then use that table to filter out outliers and create a
        // new table. We do this iteratively until no outlier is found.
        JitterCountsTable outlierTestTable = null;
        boolean outlierFound = false;
        Set<MicrosatelliteSiteAnalyser> outliers = Collections.newSetFromMap(new IdentityHashMap<>());

        while(true)
        {
            JitterCountsTable newTable = new JitterCountsTable(repeatUnit, consensusType);

            for(MicrosatelliteSiteAnalyser microsatelliteSiteAnalyser : microsatelliteSiteAnalysers)
            {
                if(outliers.contains(microsatelliteSiteAnalyser))
                {
                    continue;
                }

                if(!microsatelliteSiteAnalyser.shouldKeepSite(JitterAnalyserConstants.ALT_COUNT_FRACTION_INIT,
                        JitterAnalyserConstants.ALT_COUNT_FRACTION_STEP,
                        JitterAnalyserConstants.MAX_REJECTED_READ_FRACTION,
                        JitterAnalyserConstants.MIN_PASSING_SITE_READS))
                {
                    continue;
                }

                // get all the read counts into a row object
                JitterTableRow row = new JitterTableRow(
                        microsatelliteSiteAnalyser.refGenomeMicrosatellite().numRepeat, newTable.RepeatUnit, newTable.ConsensusType);

                for(Map.Entry<Integer, Integer> entry : microsatelliteSiteAnalyser.passingJitterCounts(consensusType).entrySet())
                {
                    int jitter = entry.getKey();
                    int numReads = entry.getValue();
                    row.addReads(jitter, numReads);
                }

                newTable.mergeCounts(row);
            }

            if(outlierTestTable != null && !outlierFound)
            {
                return newTable;
            }

            outlierFound = false;
            outlierTestTable = newTable;
        }
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

        return selectors;
    }

    private List<JitterModelParams> fitJitterModels(final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers)
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

        List<JitterModelParams> fittedParams = Lists.newArrayList();
        for(ConsensusType consensusType : consensusTypes())
        {
            for(MicrosatelliteSelector selector : selectors)
            {
                JitterCountsTable msStatsTable = buildJitterCountsTable(
                        selector.unitName(), consensusType,
                        microsatelliteSiteAnalysers.stream().filter(selector::select).collect(Collectors.toList()));

                JitterModelFitter fitter = new JitterModelFitter(msStatsTable);
                fitter.performFit();
                fittedParams.add(fitter.getJitterModelParams());
            }
        }

        return fittedParams;
    }

    private static void drawMicrosatelliteCharts(final String outputDir, final String sampleId, final String statsTableFile)
            throws IOException, InterruptedException
    {
        int result = RExecutor.executeFromClasspath("basequal/msi_jitter_plot.R", outputDir, sampleId, statsTableFile);
        if(result != 0)
            throw new IOException("R execution failed. Unable to complete segmentation.");
    }
}
