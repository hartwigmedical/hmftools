package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.targeted.TargetRegions.mergeTumorRatios;

import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.TargetRegions;
import com.hartwig.hmftools.cobalt.targeted.WholeGenome;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class CobaltCalculator
{
    private final ListMultimap<HumanChromosome, CobaltRatio> mRatios;
    private final List<MedianRatio> mMedianRatios;
    private final GcMedianReadDepth mTumorStats;
    private final GcMedianReadDepth mReferenceStatistics;

    public CobaltCalculator(
            final ListMultimap<HumanChromosome, DepthReading> tumourDepthReadings,
            final ListMultimap<HumanChromosome, DepthReading> referenceDepthReadings,
            final CobaltConfig config)
    {
        Preconditions.checkArgument(!tumourDepthReadings.isEmpty() || !referenceDepthReadings.isEmpty());
        WindowStatuses windowStatuses = new WindowStatuses(config.gcProfileData(), config.excludedRegions(), config.diploidRegions());

        ResultsConsolidator resultsConsolidator = null;
        ListMultimap<HumanChromosome,BamRatio> tumorResults = null;
        ListMultimap<HumanChromosome, BamRatio> referenceResults = null;

        if(config.targetedPanelMode())
        {
            List<TargetRegions> panelScopes = config.targetRegionScopes();

            if(panelScopes.isEmpty())
            {
                CB_LOGGER.error("invalid target region norm file scope state");
                System.exit(1);
            }

            tumorResults = ArrayListMultimap.create();
            referenceResults = ArrayListMultimap.create();

            GcMedianReadDepth tumorGcMedianReadDepth = null;
            List<MedianRatio> referenceMedianRatios = null;
            GcMedianReadDepth refGcMedianReadDepth = null;

            for(TargetRegions targetRegionScope : panelScopes)
            {
                TumorCalculation tumorCalc = new TumorCalculation(windowStatuses, targetRegionScope);
                tumourDepthReadings.forEach(tumorCalc::addReading);
                ListMultimap<HumanChromosome, BamRatio> panelRatios = tumorCalc.calculateRatios();

                mergeTumorRatios(tumorResults, panelRatios);

                if(resultsConsolidator == null)
                {
                    // only the first of these are used when multiple panels are handled
                    resultsConsolidator = tumorCalc.consolidator();

                    // only the first panel's GC read info will be written to file
                    tumorGcMedianReadDepth = tumorCalc.medianReadDepths();
                }

                ReferenceCalculation referenceCalc = new ReferenceCalculation(
                        windowStatuses, targetRegionScope, config.refGenomeVersion(), resultsConsolidator, !config.targetedPanelMode());

                referenceDepthReadings.forEach(referenceCalc::addReading);
                ListMultimap<HumanChromosome, BamRatio> refPanelRatios = referenceCalc.calculateRatios();

                mergeTumorRatios(referenceResults, refPanelRatios);

                GcMedianReadDepth refTargetRegionGcMedianReadDepth = referenceCalc.medianReadDepths();

                if(referenceMedianRatios == null)
                {
                    // as per tumor, only take the initial target region values
                    referenceMedianRatios = referenceCalc.medianRatios(referenceResults, config.refGenomeVersion());
                    refGcMedianReadDepth = refTargetRegionGcMedianReadDepth;
                }

                if(!referenceDepthReadings.isEmpty())
                {
                    CB_LOGGER.info("reference sample median({}), mean({})",
                            formatReadDepth(refTargetRegionGcMedianReadDepth.medianReadDepth()), formatReadDepth(refTargetRegionGcMedianReadDepth.meanReadDepth()));
                }
            }

            mTumorStats = tumorGcMedianReadDepth;
            mMedianRatios = referenceMedianRatios;
            mReferenceStatistics = refGcMedianReadDepth;
        }
        else
        {
            WholeGenome scope = new WholeGenome();

            TumorCalculation tumorCalc = new TumorCalculation(windowStatuses, scope);
            tumourDepthReadings.forEach(tumorCalc::addReading);
            tumorResults = tumorCalc.calculateRatios();
            mTumorStats = tumorCalc.medianReadDepths();

            resultsConsolidator = tumorCalc.consolidator();

            ReferenceCalculation referenceCalc = new ReferenceCalculation(
                    windowStatuses, scope, config.refGenomeVersion(), resultsConsolidator, !config.targetedPanelMode());

            referenceDepthReadings.forEach(referenceCalc::addReading);
            referenceResults = referenceCalc.calculateRatios();
            mMedianRatios = referenceCalc.medianRatios(referenceResults, config.refGenomeVersion());
            mReferenceStatistics = referenceCalc.medianReadDepths();

            if(!referenceDepthReadings.isEmpty())
            {
                CB_LOGGER.info("reference sample median({}), mean({})",
                        formatReadDepth(mReferenceStatistics.medianReadDepth()), formatReadDepth(mReferenceStatistics.meanReadDepth()));
            }
        }

        if(!tumourDepthReadings.isEmpty())
        {
            CB_LOGGER.info("tumor sample median({}), mean({}}",
                    formatReadDepth(mTumorStats.medianReadDepth()), formatReadDepth(mTumorStats.meanReadDepth()));
        }

        ResultsCollator collator = new ResultsCollator(config.refGenomeVersion());
        mRatios = collator.collateResults(tumorResults, referenceResults);
    }

    public ListMultimap<HumanChromosome, CobaltRatio> getCalculatedRatios() { return mRatios; }

    public List<MedianRatio> medianRatios() { return mMedianRatios; }
    public GcMedianReadDepth tumorMedianReadDepth() { return mTumorStats; }
    public GcMedianReadDepth referenceMedianReadDepth() { return mReferenceStatistics; }

    private static String formatReadDepth(Double value)
    {
        return String.format("%.2f", value);
    }
}
