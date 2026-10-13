package com.hartwig.hmftools.cobalt.calculations;

import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.targeted.TargetRegions.mergeTumorRatios;

import java.util.Collections;
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

        List<MedianRatio> referenceMedianRatios = null;
        GcMedianReadDepth refGcMedianReadDepth = null;
        ListMultimap<HumanChromosome, BamRatio> referenceResults = null;

        boolean hasReference = config.hasReferenceId();

        if(!hasReference)
        {
            // initialise with empty values
            referenceMedianRatios = Collections.emptyList();
            refGcMedianReadDepth = GcMedianReadDepth.NO_RESULTS;
            referenceResults = ArrayListMultimap.create();
        }

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

            for(int i = 0; i < panelScopes.size(); ++i)
            {
                TargetRegions targetRegionScope = panelScopes.get(i);
                TumorCalculation tumorCalc = new TumorCalculation(windowStatuses, targetRegionScope);

                tumourDepthReadings.forEach(tumorCalc::addReading);

                ListMultimap<HumanChromosome, BamRatio> panelRatios = tumorCalc.calculateRatios();

                logReadDepthInfo(tumorCalc, format("tumor panel %d", i));

                mergeTumorRatios(tumorResults, panelRatios);

                if(resultsConsolidator == null)
                    resultsConsolidator = tumorCalc.consolidator();

                if(tumorGcMedianReadDepth == null)
                {
                    // only the first panel's GC read info will be cached and written to file
                    tumorGcMedianReadDepth = tumorCalc.medianReadDepths();
                }

                if(hasReference)
                {
                    ReferenceCalculation referenceCalc = new ReferenceCalculation(
                            windowStatuses, targetRegionScope, config.refGenomeVersion(), resultsConsolidator, !config.targetedPanelMode());

                    referenceDepthReadings.forEach(referenceCalc::addReading);

                    ListMultimap<HumanChromosome, BamRatio> refPanelRatios = referenceCalc.calculateRatios();

                    logReadDepthInfo(referenceCalc, format("reference panel %d", i));

                    mergeTumorRatios(referenceResults, refPanelRatios);

                    GcMedianReadDepth refTargetRegionGcMedianReadDepth = referenceCalc.medianReadDepths();

                    if(referenceMedianRatios == null)
                    {
                        // as per tumor, only take the initial target region values
                        referenceMedianRatios = referenceCalc.medianRatios(referenceResults, config.refGenomeVersion());
                        refGcMedianReadDepth = refTargetRegionGcMedianReadDepth;
                    }
                }
            }

            mTumorStats = tumorGcMedianReadDepth;
        }
        else
        {
            WholeGenome scope = new WholeGenome();

            TumorCalculation tumorCalc = new TumorCalculation(windowStatuses, scope);

            tumourDepthReadings.forEach(tumorCalc::addReading);

            tumorResults = tumorCalc.calculateRatios();

            logReadDepthInfo(tumorCalc, "tumor");

            mTumorStats = tumorCalc.medianReadDepths();

            resultsConsolidator = tumorCalc.consolidator();

            if(hasReference)
            {
                ReferenceCalculation referenceCalc = new ReferenceCalculation(
                        windowStatuses, scope, config.refGenomeVersion(), resultsConsolidator, !config.targetedPanelMode());

                referenceDepthReadings.forEach(referenceCalc::addReading);

                referenceResults = referenceCalc.calculateRatios();

                logReadDepthInfo(referenceCalc, "reference");

                referenceMedianRatios = referenceCalc.medianRatios(referenceResults, config.refGenomeVersion());
                refGcMedianReadDepth = referenceCalc.medianReadDepths();
            }
        }

        mMedianRatios = referenceMedianRatios;
        mReferenceStatistics = refGcMedianReadDepth;

        ResultsCollator collator = new ResultsCollator(config.refGenomeVersion());
        mRatios = collator.collateResults(tumorResults, referenceResults);
    }

    public ListMultimap<HumanChromosome, CobaltRatio> getCalculatedRatios() { return mRatios; }
    public List<MedianRatio> medianRatios() { return mMedianRatios; }
    public GcMedianReadDepth tumorMedianReadDepth() { return mTumorStats; }
    public GcMedianReadDepth referenceMedianReadDepth() { return mReferenceStatistics; }

    private void logReadDepthInfo(final BamCalculation bamCalculation, final String id)
    {
        GcMedianReadDepth gcMedianReadDepth = bamCalculation.medianReadDepths();

        CB_LOGGER.info(format("%s sample median(%.2f) mean(%.2f)",
                id, gcMedianReadDepth.medianReadDepth(), gcMedianReadDepth.meanReadDepth()));
    }
}
