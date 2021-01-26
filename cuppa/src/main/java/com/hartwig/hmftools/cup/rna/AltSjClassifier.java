package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.cup.common.CupCalcs.adjustLowProbabilities;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.ClassifierType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

public class AltSjClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final Map<String,List<AltSjPrevData>> mSampleAltSJs; // currently unused
    private final Map<String,List<AltSjPrevData>> mRefAltSjPrevalence; // ref alt-SJs by chromosome
    private final Set<Integer> mRefAltSjStartPositions; // for fast sample data filtering
    private final SampleDataCache mSampleDataCache;
    private boolean mIsValid;

    private static final double ZERO_PREVALENCE_ALLOCATION = 0.03;

    public AltSjClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleAltSJs = Maps.newHashMap();
        mRefAltSjPrevalence = Maps.newHashMap();
        mSampleDataCache = sampleDataCache;
        mRefAltSjStartPositions = Sets.newHashSet();
        mIsValid = true;

        if(config.RefAltSjPrevFile.isEmpty())
            return;

        mIsValid &= loadRefAltSJs(config.RefAltSjPrevFile);
    }

    public CategoryType categoryType() { return ALT_SJ; }
    public boolean isValid() { return mIsValid; }

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(!mIsValid || mRefAltSjPrevalence.isEmpty())
            return;

        final List<AltSjPrevData> sampleAltSJs = loadSampleAltSJs(sample.Id);

        if(sampleAltSJs == null || sampleAltSJs.isEmpty())
            return;

        calcCancerTypeProbability(sample, sampleAltSJs, results);
    }

    private AltSjPrevData findRefAltSjData(final AltSjPrevData altSJ)
    {
        final List<AltSjPrevData> refAltSJs = mRefAltSjPrevalence.get(altSJ.Location.Chromosome);

        if(refAltSJs == null)
            return null;

        return refAltSJs.stream().filter(x -> x.matches(altSJ)).findFirst().orElse(null);
    }

    private void calcCancerTypeProbability(
            final SampleData sample, final List<AltSjPrevData> sampleAltSJs, final List<SampleResult> results)
    {
        // taking the set of drivers as a group, report on the combined probability for each cancer type
        final Map<String,Double> cancerProbTotals = Maps.newHashMap();

        final Set<String> cancerTypes = mSampleDataCache.RefCancerSampleData.keySet();
        int matchedAltSJs = 0;

        for(final AltSjPrevData altSJ : sampleAltSJs)
        {
            // check the alt-SJ is one of the reference ones
            final AltSjPrevData refAltSJ = findRefAltSjData(altSJ);

            if(refAltSJ == null)
                continue;

            ++matchedAltSJs;

            for(String cancerType : cancerTypes)
            {
                if(!checkIsValidCancerType(sample, cancerType, cancerProbTotals))
                    continue;

                Double prevalence = refAltSJ.CancerPrevalences.get(cancerType);

                Double probabilityTotal = cancerProbTotals.get(cancerType);

                if(probabilityTotal == null)
                    probabilityTotal = 1.0;

                boolean adjustMatchingCancerPrev = sample.CancerType.equals(cancerType);

                double driverPrevValue;
                double prevalenceTotal = refAltSJ.PrevalenceTotal;

                if(prevalence != null)
                {
                    driverPrevValue = prevalence;

                    if(adjustMatchingCancerPrev)
                    {
                        int cohortSize = mSampleDataCache.getCancerSampleCount(cancerType);
                        double adjustedIncidence = max(driverPrevValue * cohortSize - 1, 0.0);
                        double adjustedDriverPrevValue = cohortSize > 1 ? adjustedIncidence / (cohortSize - 1) : 0;
                        prevalenceTotal -= driverPrevValue - adjustedDriverPrevValue;
                        driverPrevValue = adjustedDriverPrevValue;
                    }
                }
                else
                {
                    driverPrevValue = refAltSJ.MinPrevalence;
                }

                probabilityTotal *= driverPrevValue / prevalenceTotal;
                cancerProbTotals.put(cancerType, probabilityTotal);
            }

            adjustLowProbabilities(cancerProbTotals);
        }

        if(matchedAltSJs > 0)
        {
            final String altSjStr = String.format("ALT-SJs=%d", matchedAltSJs);

            SampleResult result = new SampleResult(
                    sample.Id, ALT_SJ, LIKELIHOOD, ClassifierType.ALT_SJ.toString(), altSjStr, cancerProbTotals);

            results.add(result);
        }
    }

    private boolean loadRefAltSJs(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), ",");
            lines.remove(0);

            // CancerType,GeneId,Chromosome,SjStart,SjEnd,Type,Prev
            int cancerIndex = fieldsIndexMap.get("CancerType");
            int geneIdIndex = fieldsIndexMap.get("GeneId");
            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int posStartIndex = fieldsIndexMap.get("SjStart");
            int posEndIndex = fieldsIndexMap.get("SjEnd");
            Integer typeIndex = fieldsIndexMap.get("Type"); // for info sake only
            int prevIndex = fieldsIndexMap.get("Prev");

            double noPrevalence = ZERO_PREVALENCE_ALLOCATION / mSampleDataCache.RefCancerSampleData.size();

            for(String data : lines)
            {
                final String items[] = data.split(",", -1);

                String cancerType = items[cancerIndex];
                String chromosome = items[chromosomeIndex];
                String asjType = typeIndex != null ? items[typeIndex] : "N/A";
                int asjPosStart = Integer.parseInt(items[posStartIndex]);
                int asjPosEnd = Integer.parseInt(items[posEndIndex]);
                double prevalence = Double.parseDouble(items[prevIndex]);

                AltSjPrevData altSJ = new AltSjPrevData(items[geneIdIndex], asjType, chromosome, asjPosStart, asjPosEnd, 0);

                AltSjPrevData matchedAltSJ = null;

                List<AltSjPrevData> altSJs = mRefAltSjPrevalence.get(chromosome);

                if(altSJs == null)
                {
                    altSJs = Lists.newArrayList();
                    mRefAltSjPrevalence.put(chromosome, altSJs);
                }
                else
                {
                    matchedAltSJ = altSJs.stream().filter(x -> x.matches(altSJ)).findFirst().orElse(null);
                }

                if(matchedAltSJ == null)
                {
                    altSJs.add(altSJ);
                    matchedAltSJ = altSJ;
                }

                mRefAltSjStartPositions.add(asjPosStart);

                double adjPrevalence = prevalence + noPrevalence;
                matchedAltSJ.CancerPrevalences.put(cancerType, adjPrevalence);
                matchedAltSJ.PrevalenceTotal += adjPrevalence;
                matchedAltSJ.MinPrevalence = noPrevalence;
            }

            // set default value for missing cancer types
            final Set<String> cancerTypes = mSampleDataCache.RefCancerSampleData.keySet();

            for(List<AltSjPrevData> altSJs : mRefAltSjPrevalence.values())
            {
                for(AltSjPrevData altSJ : altSJs)
                {
                    for(String cancerType : cancerTypes)
                    {
                        if(altSJ.CancerPrevalences.containsKey(cancerType))
                            continue;

                        altSJ.PrevalenceTotal += noPrevalence;
                    }
                }
            }

            CUP_LOGGER.info("loaded {} ref alt splice junctions", mRefAltSjPrevalence.values().stream().mapToInt(x -> x.size()).sum());
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to load ref alt splice junction file({}): {}", filename.toString(), e.toString());
            return false;
        }

        return true;
    }

    private static final String ALT_SJ_FILE_ID = ".isf.alt_splice_junc.csv";
    private static final int ALT_SJ_FRAG_COUNT_THRESHOLD = 3;

    private List<AltSjPrevData> loadSampleAltSJs(final String sampleId)
    {
        final String filename = mConfig.SampleDataDir + sampleId + ALT_SJ_FILE_ID;

        final List<AltSjPrevData> sampleAltSJs = Lists.newArrayList();

        if(!Files.exists(Paths.get(filename)))
            return sampleAltSJs;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), ",");
            lines.remove(0);

            int geneIdIndex = fieldsIndexMap.get("GeneId");
            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int posStartIndex = fieldsIndexMap.get("SjStart");
            int posEndIndex = fieldsIndexMap.get("SjEnd");
            int typeIndex = fieldsIndexMap.get("Type");
            int fragCountIndex = fieldsIndexMap.get("FragCount");

            for(String data : lines)
            {
                final String items[] = data.split(",", -1);

                int asjPosStart = Integer.parseInt(items[posStartIndex]);

                if(!mRefAltSjStartPositions.contains(asjPosStart))
                    continue;

                int fragCount = Integer.parseInt(items[fragCountIndex]);

                if(fragCount < ALT_SJ_FRAG_COUNT_THRESHOLD)
                    continue;

                String chromosome = items[chromosomeIndex];
                int asjPosEnd = Integer.parseInt(items[posEndIndex]);

                sampleAltSJs.add(new AltSjPrevData(items[geneIdIndex], items[typeIndex], chromosome, asjPosStart, asjPosEnd));
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to load alt splice junction file({}): {}", filename.toString(), e.toString());
            return null;
        }

        return sampleAltSJs;
    }
}
