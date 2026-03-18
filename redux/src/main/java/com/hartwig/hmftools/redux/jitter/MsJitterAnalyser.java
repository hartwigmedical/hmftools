package com.hartwig.hmftools.redux.jitter;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.DUAL_BASE_1;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.DUAL_BASE_2;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.DUAL_BASE_3;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.DUAL_BASE_4;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.SINGLE_BASE_1;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.SINGLE_BASE_2;
import static com.hartwig.hmftools.redux.ms_model.MsModelParams.DEFAULT_MODEL_PARAMS;

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
import com.hartwig.hmftools.common.redux.JitterCountsTable;
import com.hartwig.hmftools.common.redux.JitterCountsTableFile;
import com.hartwig.hmftools.common.redux.JitterModelParams;
import com.hartwig.hmftools.common.redux.JitterModelParamsFile;
import com.hartwig.hmftools.common.redux.JitterTableRow;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.redux.MsiModelPrediction;
import com.hartwig.hmftools.common.utils.RExecutor;
import com.hartwig.hmftools.redux.ms_model.MsModelCalculator;

import htsjdk.samtools.SAMRecord;

public class MsJitterAnalyser
{
    private final MsJitterConfig mConfig;
    private final SampleReadProcessor mSampleReadProcessor;

    private EnumSet<ConsensusType> mConsensusTypes;

    public MsJitterAnalyser(final MsJitterConfig config)
    {
        mConfig = config;

        List<MicrosatelliteSite> microsatelliteSites = loadRefGenomeMicrosatellites();
        ConsensusMarker consensusMarker = ConsensusMarker.create();
        mSampleReadProcessor = new SampleReadProcessor(config, microsatelliteSites, consensusMarker);

        mConsensusTypes = null;
    }

    public void processRead(final SAMRecord read)
    {
        mSampleReadProcessor.processRead(read);
    }

    private EnumSet<ConsensusType> consensusTypes()
    {
        if(mConsensusTypes != null)
            return mConsensusTypes;

        mConsensusTypes = MsJitterConfig.consensusTypes(mConfig);
        return mConsensusTypes;
    }

    public void writeAnalysisOutput() throws IOException, InterruptedException
    {
        Collection<MicrosatelliteSiteData> microsatelliteSiteData = mSampleReadProcessor.getMicrosatelliteSiteAnalysers();

        // now write out all the repeat stats
        if(mConfig.WriteSiteFile)
        {
            MicrosatelliteSiteFile.write(MicrosatelliteSiteFile.generateFilename(
                    mConfig.OutputDir, mConfig.SampleId), microsatelliteSiteData, consensusTypes());
        }

        List<JitterCountsTable> jitterCountsTable = Lists.newArrayList();

        for(ConsensusType consensusType : consensusTypes())
        {
            for(MicrosatelliteSelector s : createMicrosatelliteSelectorsForCharts())
            {
                JitterCountsTable msStatsTable = buildJitterCountsTable(
                        s.unitName(), consensusType, microsatelliteSiteData.stream().filter(s::select).collect(Collectors.toList()));

                jitterCountsTable.add(msStatsTable);
            }
        }

        writeMicrosatelliteStatsTable(jitterCountsTable);

        if(mConfig.runMsiPrediction())
            writeMsModelPrediction(jitterCountsTable);

        // now perform the fitting
        List<JitterModelParams> jitterModelParamsList = fitJitterModels(microsatelliteSiteData);

        JitterModelParamsFile.write(JitterModelParamsFile.generateFilename(mConfig.OutputDir, mConfig.SampleId), jitterModelParamsList);
    }

    private void writeMsModelPrediction(final List<JitterCountsTable> jitterCountsTable) throws IOException
    {
        MsModelCalculator modelCalculator = new MsModelCalculator(
                DEFAULT_MODEL_PARAMS, mConfig.MsModelCoefficientsFile, mConfig.MsModelErroRatesFile);

        double predictedValue = modelCalculator.calcPredictedMsiRate(jitterCountsTable);

        String filename = MsiModelPrediction.generateFilename(mConfig.OutputDir, mConfig.SampleId);
        MsiModelPrediction.write(filename, predictedValue);
    }

    private List<MicrosatelliteSite> loadRefGenomeMicrosatellites()
    {
        List<MicrosatelliteSite> microsatelliteSites = RefGenomeMicrosatelliteFile.read(mConfig.RefGenomeMsiFile);

        if(mConfig.SpecificChrRegions.hasFilters())
        {
            microsatelliteSites = microsatelliteSites.stream()
                    .filter(x -> mConfig.SpecificChrRegions.includeRegion(x.Region)).collect(Collectors.toList());
        }

        RD_LOGGER.info("loaded {} microsatellites regions", microsatelliteSites.size());
        return microsatelliteSites;
    }

    private void writeMicrosatelliteStatsTable(List<JitterCountsTable> jitterCountsTable) throws IOException, InterruptedException
    {
        String statsTableFile = JitterCountsTableFile.generateFilename(mConfig.OutputDir, mConfig.SampleId);

        JitterCountsTableFile.write(statsTableFile, jitterCountsTable);

        // draw a chart of the 9 ms profiles
        if(mConfig.WritePlots)
            drawMicrosatelliteCharts(statsTableFile);
    }

    static JitterCountsTable buildJitterCountsTable(
            final String repeatUnit, final ConsensusType consensusType,
            final Collection<MicrosatelliteSiteData> microsatelliteSites)
    {
        // In order to filter out outliers, we perform the stats summation in a loop
        // We create a table of all read stats, then use that table to filter out outliers and create a
        // new table. We do this iteratively until no outlier is found.
        JitterCountsTable outlierTestTable = null;
        boolean outlierFound = false;
        Set<MicrosatelliteSiteData> outliers = Collections.newSetFromMap(new IdentityHashMap<>());

        while(true)
        {
            JitterCountsTable newTable = new JitterCountsTable(repeatUnit, consensusType);

            for(MicrosatelliteSiteData microsatelliteSiteData : microsatelliteSites)
            {
                if(outliers.contains(microsatelliteSiteData))
                {
                    continue;
                }

                if(!microsatelliteSiteData.shouldKeepSite())
                {
                    continue;
                }

                // get all the read counts into a row object
                JitterTableRow row = new JitterTableRow(
                        microsatelliteSiteData.refGenomeMicrosatellite().RepeatCount, newTable.RepeatUnit, newTable.ConsensusType);

                for(Map.Entry<Integer, Integer> entry : microsatelliteSiteData.passingJitterCounts(consensusType).entrySet())
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

    private List<JitterModelParams> fitJitterModels(final Collection<MicrosatelliteSiteData> microsatelliteSiteData)
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
                        microsatelliteSiteData.stream().filter(selector::select).collect(Collectors.toList()));

                JitterModelFitter fitter = new JitterModelFitter(msStatsTable);
                fitter.performFit();
                fittedParams.add(fitter.getJitterModelParams());
            }
        }

        return fittedParams;
    }

    private void drawMicrosatelliteCharts(final String statsTableFile) throws IOException, InterruptedException
    {
        int result = RExecutor.executeFromClasspath("basequal/msi_jitter_plot.R", mConfig.OutputDir, mConfig.SampleId, statsTableFile);
        if(result != 0)
            throw new IOException("R execution failed. Unable to complete segmentation.");
    }
}
