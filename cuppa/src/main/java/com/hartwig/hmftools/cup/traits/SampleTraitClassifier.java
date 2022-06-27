package com.hartwig.hmftools.cup.traits;

import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcPercentilePrevalence;
import static com.hartwig.hmftools.cup.common.CupConstants.UNDEFINED_PERC_MAX_MULTIPLE;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PERCENTILE;
import static com.hartwig.hmftools.cup.common.ResultType.PREVALENCE;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.GENDER;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.MS_INDELS_TMB;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.PLOIDY;
import static com.hartwig.hmftools.cup.traits.SampleTraitType.WGD;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.GENDER_FEMALE_INDEX;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.GENDER_MALE_INDEX;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.loadRefGenderRateData;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.loadTraitsFromCohortFile;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.loadTraitsFromDatabase;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.loadRefPercentileData;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.loadRefRateData;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.ClassifierType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.ResultType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SampleTraitClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,SampleTraitsData> mSampleTraitsData;

    private final Map<SampleTraitType,Map<String,double[]>> mRefTraitPercentiles;
    private final Map<SampleTraitType,Map<String,Double>> mRefTraitRates;
    private final boolean mApplyPloidyLikelihood;

    private final Map<String,double[]> mRefGenderRates;

    private static final String APPLY_PLOIDY_LIKELIHOOD = "apply_ploidy_likelihood";

    public SampleTraitClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleTraitsData = Maps.newHashMap();
        mRefTraitPercentiles = Maps.newHashMap();
        mRefTraitRates = Maps.newHashMap();
        mRefGenderRates = Maps.newHashMap();

        mApplyPloidyLikelihood = cmd.hasOption(APPLY_PLOIDY_LIKELIHOOD);
    }

    public CategoryType categoryType() { return SAMPLE_TRAIT; }
    public void close() {}

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(APPLY_PLOIDY_LIKELIHOOD, false, "Add ploidy high/low likelihood feature");
    }

    @Override
    public boolean loadData()
    {
        return loadRefData() && loadSampleData();
    }

    private boolean loadRefData()
    {
        if(mConfig.RefTraitPercFile.isEmpty() && mConfig.RefTraitRateFile.isEmpty())
            return false;

        if(!loadRefPercentileData(mConfig.RefTraitPercFile, mRefTraitPercentiles))
        {
            CUP_LOGGER.error("invalid ref sample trait percentile data");
            return false;
        }

        if(!loadRefRateData(mConfig.RefTraitRateFile, mRefTraitRates))
        {
            CUP_LOGGER.error("invalid ref sample trait rates data");
            return false;
        }

        if(Files.exists(Paths.get(mConfig.RefGenderRateFile)) && !loadRefGenderRateData(mConfig.RefGenderRateFile, mRefGenderRates))
        {
            CUP_LOGGER.error("invalid ref gender rates data");
            return false;
        }

        return true;
    }

    private boolean loadSampleData()
    {
        if(mConfig.TestRefData)
        {
            if(!mConfig.RefSampleTraitsFile.isEmpty())
            {
                CUP_LOGGER.info("loading ref cohort traits from file({})", mConfig.RefSampleTraitsFile);

                if(!loadTraitsFromCohortFile(mConfig.RefSampleTraitsFile, mSampleTraitsData))
                    return false;

                CUP_LOGGER.info("loaded traits for {} samples", mSampleTraitsData.size());
            }
            else
            {
                CUP_LOGGER.error("missing ref cohort traits file");
                return false;
            }
        }
        else
        {
            if(mConfig.DbAccess != null)
            {
                if(!loadTraitsFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleTraitsData))
                    return false;
            }
            else
            {
                for(SampleData sample : mSampleDataCache.SampleDataList)
                {
                    String purpleDir = mConfig.getPurpleDataDir(sample.Id);

                    try
                    {
                        final PurityContext purityContext = PurityContextFile.read(purpleDir, sample.Id);
                        SampleTraitsData traitsData = SampleTraitsData.from(sample.Id, purityContext, 0);
                        mSampleTraitsData.put(traitsData.SampleId, traitsData);
                    }
                    catch(Exception e)
                    {
                        CUP_LOGGER.error("sample({}) sample traits - failed to load purity file from dir{}): {}",
                                sample.Id, purpleDir, e.toString());
                        return false;
                    }
                }
            }
        }

        // cache traits for other components to use
        for(SampleData sample : mSampleDataCache.SampleDataList)
        {
            final SampleTraitsData sampleTraits = mSampleTraitsData.get(sample.Id);

            if(sampleTraits != null)
                sample.setGender(sampleTraits.GenderType);
        }

        mSampleTraitsData.entrySet().forEach(x -> mSampleDataCache.SampleTraitsData.put(x.getKey(), x.getValue()));

        return true;
    }

    @Override
    public boolean processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(mRefTraitRates.isEmpty())
            return false;

        final SampleTraitsData sampleTraits = mSampleTraitsData.get(sample.Id);

        if(sampleTraits == null)
        {
            CUP_LOGGER.warn("sample({}) has missing traits data, classifier invalid", sample.Id);
            return false;
        }

        addTraitPrevalences(sample, sampleTraits, results);

        int cancerSampleCount = sample.isRefSample() ? mSampleDataCache.getCancerSampleCount(sample.cancerType()) : 0;
        addTraitLikelihoods(sample, cancerSampleCount, results, MS_INDELS_TMB, sampleTraits.IndelsMbPerMb);

        if(mApplyPloidyLikelihood)
            addTraitLikelihoods(sample, cancerSampleCount, results, PLOIDY, sampleTraits.Ploidy);

        addGenderClassifier(sample, sampleTraits, results);
        return true;
    }

    private void addGenderClassifier(final SampleData sample, final SampleTraitsData sampleTraits, final List<SampleResult> results)
    {
        if(mRefGenderRates.isEmpty())
            return;

        final Map<String, Double> cancerProbs = Maps.newHashMap();

        for(Map.Entry<String,double[]> entry : mRefGenderRates.entrySet())
        {
            final String cancerType = entry.getKey();
            final double[] genderRates = entry.getValue();

            if(sampleTraits.GenderType == Gender.FEMALE)
                cancerProbs.put(cancerType, genderRates[GENDER_FEMALE_INDEX]);
            else
                cancerProbs.put(cancerType, genderRates[GENDER_MALE_INDEX]);
        }

        SampleResult result = new SampleResult(
                sample.Id, SAMPLE_TRAIT, ResultType.CLASSIFIER, ClassifierType.GENDER.toString(), sampleTraits.GenderType.toString(), cancerProbs);

        results.add(result);
    }

    private void addTraitPrevalences(final SampleData sample, final SampleTraitsData sampleTraits, final List<SampleResult> results)
    {
        for(Map.Entry<SampleTraitType,Map<String,Double>> entry : mRefTraitRates.entrySet())
        {
            final SampleTraitType traitType = entry.getKey();

            if(!isReportableType(traitType))
                continue;

            Map<String, Double> cancerRates = entry.getValue();

            // reverse the prevalence for MALE since gender is currently IsFemale
            if(traitType == GENDER)
            {
                if(sampleTraits.GenderType != Gender.FEMALE)
                {
                    Map<String, Double> oppGenderRates = Maps.newHashMap();
                    cancerRates.entrySet().forEach(x -> oppGenderRates.put(x.getKey(), 1 - x.getValue()));
                    cancerRates = oppGenderRates;
                }

                SampleResult result = new SampleResult(
                        sample.Id, SAMPLE_TRAIT, PREVALENCE, traitType.toString(), sampleTraits.getStrValue(traitType), cancerRates);

                results.add(result);
            }
            else if(traitType == WGD)
            {
                SampleResult result = new SampleResult(
                        sample.Id, SAMPLE_TRAIT, PREVALENCE, traitType.toString(), sampleTraits.getStrValue(traitType), cancerRates);

                results.add(result);
            }
        }

        for(Map.Entry<SampleTraitType, Map<String, double[]>> entry : mRefTraitPercentiles.entrySet())
        {
            final SampleTraitType traitType = entry.getKey();

            if(!isReportableType(traitType))
                continue;

            double traitValue = sampleTraits.getDoubleValue(traitType);

            final Map<String, Double> cancerTypeValues = Maps.newHashMap();

            for(Map.Entry<String, double[]> cancerPercentiles : entry.getValue().entrySet())
            {
                final String cancerType = cancerPercentiles.getKey();

                if(!isKnownCancerType(cancerType))
                    continue;

                double percentile = getPercentile(cancerPercentiles.getValue(), traitValue, true, UNDEFINED_PERC_MAX_MULTIPLE);
                cancerTypeValues.put(cancerType, percentile);
            }

            SampleResult result = new SampleResult(
                    sample.Id, SAMPLE_TRAIT, PERCENTILE, traitType.toString(), String.valueOf(traitValue), cancerTypeValues);

            results.add(result);
        }
    }

    private void addTraitLikelihoods(
            final SampleData sample, int cancerSampleCount, final List<SampleResult> results, final SampleTraitType traitType, double traitValue)
    {
        int cancerTypeCount = mSampleDataCache.RefCancerSampleData.size();

        final Map<String, double[]> indelPercentiles = mRefTraitPercentiles.get(traitType);

        final Map<String, Double> cancerPrevsLow = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, indelPercentiles, traitValue, true);
        results.add(new SampleResult(sample.Id, SAMPLE_TRAIT, LIKELIHOOD, traitType.toString() + "_LOW",
                String.format("%.4f", traitValue), cancerPrevsLow));

        final Map<String, Double> cancerPrevsHigh = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, indelPercentiles, traitValue, false);
        results.add(new SampleResult(sample.Id, SAMPLE_TRAIT, LIKELIHOOD, traitType.toString() + "_HIGH",
                String.format("%.4f", traitValue), cancerPrevsHigh));
    }

    private static boolean isReportableType(final SampleTraitType type)
    {
        return (type == MS_INDELS_TMB || type == GENDER);
    }

}
