package com.hartwig.hmftools.cup.traits;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.stats.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_CANCER_TYPE;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.SUBSET_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_TRAITS_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENDER_RATES;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_TRAIT_PERC;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_TRAIT_RATES;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CupConstants.BREAST_MALE_GENDER_RATE;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_BREAST;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_BREAST_TRIPLE_NEGATIVE;
import static com.hartwig.hmftools.cup.common.CupConstants.isCandidateCancerType;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.GENDER_RATES;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.parseFileSet;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.FLD_GENDER_FEMALE;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.FLD_GENDER_MALE;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.GENDER_FEMALE_INDEX;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.GENDER_MALE_INDEX;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.loadTraitsFromCohortFile;
import static com.hartwig.hmftools.cup.traits.SampleTraitsDataLoader.loadTraitsFromDatabase;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;
import com.hartwig.hmftools.cup.ref.RefClassifier;

import org.apache.commons.cli.CommandLine;

public class RefSampleTraits implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,List<SampleTraitsData>> mCancerTraitsData;
    private final Map<String,double[]> mGenderRates;

    private BufferedWriter mPercentilesWriter;
    private BufferedWriter mRatesWriter;

    public RefSampleTraits(final RefDataConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerTraitsData = Maps.newHashMap();
        mPercentilesWriter = null;
        mRatesWriter = null;

        mGenderRates = Maps.newHashMap();

        if(cmd.hasOption(GENDER_RATES))
        {
            String[] genderEntries = cmd.getOptionValue(GENDER_RATES).split(DATA_DELIM);
            for(String genderEntry : genderEntries)
            {
                String[] genderData = genderEntry.split(SUBSET_DELIM);

                if(genderData.length == 3)
                {
                    double[] rates = new double[2];
                    rates[GENDER_MALE_INDEX] = Double.parseDouble(genderData[1 + GENDER_MALE_INDEX]);
                    rates[GENDER_FEMALE_INDEX] = Double.parseDouble(genderData[1 + GENDER_FEMALE_INDEX]);
                    mGenderRates.put(genderData[0], rates);
                }
            }
        }
        else
        {
            // -gender_rates "Breast;1;0.1"
            mGenderRates.put(CANCER_TYPE_BREAST, new double[] {1.0, BREAST_MALE_GENDER_RATE} );
            mGenderRates.put(CANCER_TYPE_BREAST_TRIPLE_NEGATIVE, new double[] {1.0, BREAST_MALE_GENDER_RATE} );
        }

        // add in the default zero-prevalence ones
        for(String cancerType : mSampleDataCache.RefCancerSampleData.keySet())
        {
            if(mGenderRates.containsKey(cancerType))
                continue;

            double[] rates = new double[2];

            rates[GENDER_FEMALE_INDEX] = isCandidateCancerType(Gender.FEMALE, cancerType) ? 1 : 0;
            rates[GENDER_MALE_INDEX] = isCandidateCancerType(Gender.MALE, cancerType) ? 1 : 0;
            mGenderRates.put(cancerType, rates);
        }
    }

    public CategoryType categoryType() { return SAMPLE_TRAIT; }
    public static boolean requiresBuild(final RefDataConfig config)
    {
        return config.Categories.contains(SAMPLE_TRAIT) || !config.CohortSampleTraitsFile.isEmpty() || config.DbAccess != null;
    }

    public void buildRefDataSets()
    {
        if(mConfig.CohortSampleTraitsFile.isEmpty() && mConfig.DbAccess == null && mConfig.PurpleDir.isEmpty())
            return;

        CUP_LOGGER.info("building sample traits reference data");

        final Map<String,SampleTraitsData> sampleTraitsData = Maps.newHashMap();

        if(!mConfig.CohortSampleTraitsFile.isEmpty())
        {
            final List<String> files = parseFileSet(mConfig.CohortSampleTraitsFile);
            files.forEach(x -> loadRefPurityData(x, sampleTraitsData));
        }
        else if(mConfig.DbAccess != null)
        {
            loadTraitsFromDatabase(mConfig.DbAccess, mSampleDataCache.refSampleIds(false), sampleTraitsData);
            sampleTraitsData.values().forEach(x -> assignSampleTraitsData(x));
        }
        else
        {
            // load from per-sample files
            for(SampleData sample : mSampleDataCache.RefSampleDataList)
            {
                final String purpleDataDir = formSamplePath(mConfig.PurpleDir, sample.Id);

                try
                {
                    final PurityContext purityContext = PurityContextFile.read(purpleDataDir, sample.Id);

                    // CUP_LOGGER.debug("sample({}) loading sample traits from purpleDir({})", sample.Id, purpleDataDir);

                    SampleTraitsData traitsData = SampleTraitsData.from(sample.Id, purityContext, 0);
                    assignSampleTraitsData(traitsData);
                    sampleTraitsData.put(sample.Id, traitsData);
                }
                catch(Exception e)
                {
                    CUP_LOGGER.error("sample({}) sample traits - failed to load purity file from dir{}): {}",
                            sample.Id, purpleDataDir, e.toString());
                    break;
                }
            }
        }

        writeCohortData(sampleTraitsData);

        initialiseRefDataWriters();

        for(Map.Entry<String,List<SampleTraitsData>> entry : mCancerTraitsData.entrySet())
        {
            final String cancerType = entry.getKey();

            if(!isKnownCancerType(cancerType))
                continue;

            final List<SampleTraitsData> traitsData = entry.getValue();

            final List<Double> purityValues = traitsData.stream().map(x -> x.Purity).collect(Collectors.toList());
            writePercentilesData(cancerType, SampleTraitType.PURITY, createPercentileData(purityValues));

            final List<Double> ploidyValues = traitsData.stream().map(x -> x.Ploidy).collect(Collectors.toList());
            writePercentilesData(cancerType, SampleTraitType.PLOIDY, createPercentileData(ploidyValues));

            final List<Double> msIndelTmbValues = traitsData.stream().map(x -> x.IndelsMbPerMb).collect(Collectors.toList());
            writePercentilesData(cancerType, SampleTraitType.MS_INDELS_TMB, createPercentileData(msIndelTmbValues));

            double cancerSamples = mSampleDataCache.getCancerSampleCount(cancerType);
            int wgdCount = (int)traitsData.stream().filter(x -> x.HasWGD).count();
            int femaleCount = (int)traitsData.stream().filter(x -> x.GenderType == Gender.FEMALE).count();

            writeRatesData(cancerType, wgdCount/cancerSamples, femaleCount/cancerSamples);
        }

        // create a gender rates file
        writeGenderRates();

        closeBufferedWriter(mPercentilesWriter);
        closeBufferedWriter(mRatesWriter);
    }
    
    private void initialiseRefDataWriters()
    {
        try
        {
            final String percFilename = mConfig.OutputDir + REF_FILE_TRAIT_PERC;
            mPercentilesWriter = createBufferedWriter(percFilename, false);

            mPercentilesWriter.write("CancerType,TraitType");

            for(int i = 0; i < PERCENTILE_COUNT; ++i)
            {
                mPercentilesWriter.write(String.format(",Pct_%.2f", i * 0.01));
            }

            mPercentilesWriter.newLine();

            final String ratesFilename = mConfig.OutputDir + REF_FILE_TRAIT_RATES;
            mRatesWriter = createBufferedWriter(ratesFilename, false);

            mRatesWriter.write("CancerType,WGDPerc,GenderFemalePerc");
            mRatesWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample traits ref data output: {}", e.toString());
        }
    }

    private void writeRatesData(final String cancerType, double wgdRate, double femaleRate)
    {
        try
        {
            mRatesWriter.write(String.format("%s,%.4f,%.4f", cancerType, wgdRate, femaleRate));
            mRatesWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref sample traits rates ref data output: {}", e.toString());
        }
    }

    private void writePercentilesData(final String cancerType, final SampleTraitType traitType, final double[] percentileValues)
    {
        try
        {
            mPercentilesWriter.write(String.format("%s,%s", cancerType, traitType));

            for(int i = 0; i < percentileValues.length; ++i)
            {
                mPercentilesWriter.write(String.format(",%.6f", percentileValues[i]));
            }

            mPercentilesWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample traits ref data output: {}", e.toString());
        }
    }

    private double[] createPercentileData(final List<Double> values)
    {
        // sort the data into an array
        final List<Integer> sortedIndices = getSortedVectorIndices(convertList(values), true);
        final double[] sortedValues = new double[values.size()];
        
        for(int i = 0; i < values.size(); ++i)
        {
            sortedValues[i] = values.get(sortedIndices.get(i));
        }
        
        return buildPercentiles(sortedValues);
    }

    private void assignSampleTraitsData(final SampleTraitsData traitsData)
    {
        final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(traitsData.SampleId);
        if(cancerType == null)
        {
            CUP_LOGGER.error("sample({}) traits missing cancer type", traitsData.SampleId);
            return;
        }

        // cache for other components to use
        mSampleDataCache.RefSampleTraitsData.put(traitsData.SampleId, traitsData);

        if(isKnownCancerType(cancerType))
        {
            List<SampleTraitsData> traitsList = mCancerTraitsData.get(cancerType);
            if(traitsList == null)
            {
                mCancerTraitsData.put(cancerType, Lists.newArrayList(traitsData));
            }
            else
            {
                traitsList.add(traitsData);
            }
        }
    }

    private void writeGenderRates()
    {
        try
        {
            final String filename = mConfig.OutputDir + REF_FILE_GENDER_RATES;
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(String.format("%s,%s,%s", FLD_CANCER_TYPE, FLD_GENDER_FEMALE, FLD_GENDER_MALE));
            writer.newLine();

            double maleRatesTotal = mGenderRates.values().stream().mapToDouble(x -> x[GENDER_MALE_INDEX]).sum();
            double femaleRatesTotal = mGenderRates.values().stream().mapToDouble(x -> x[GENDER_FEMALE_INDEX]).sum();

            for(Map.Entry<String,double[]> entry : mGenderRates.entrySet())
            {
                final String cancerType = entry.getKey();
                final double[] rates = entry.getValue();

                writer.write(String.format("%s,%.6f,%.6f",
                        cancerType, rates[GENDER_FEMALE_INDEX] / femaleRatesTotal, rates[GENDER_MALE_INDEX] / maleRatesTotal));
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample traits cohort data output: {}", e.toString());
        }

    }

    private void writeCohortData(final Map<String,SampleTraitsData> sampleTraitsData)
    {
        if(!mConfig.WriteCohortFiles)
            return;

        final String filename = mConfig.OutputDir + COHORT_REF_TRAITS_DATA_FILE;
        if(Files.exists(Paths.get(filename)))
        {
            CUP_LOGGER.warn("not over-writing cohort sample traits reference file({})", filename);
            return;
        }

        CUP_LOGGER.info("writing cohort sample traits reference data");

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(SampleTraitsData.header());
            writer.newLine();

            for(Map.Entry<String,SampleTraitsData> entry : sampleTraitsData.entrySet())
            {
                final String sampleId = entry.getKey();
                final SampleTraitsData traitsData = entry.getValue();
                writer.write(String.format("%s,%s", sampleId, traitsData.toCsv()));
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample traits cohort data output: {}", e.toString());
        }
    }

    private void loadRefPurityData(final String filename, final Map<String,SampleTraitsData> sampleTraitsData)
    {
        final Map<String,SampleTraitsData> traitsDataMap = Maps.newHashMap();

        loadTraitsFromCohortFile(filename, traitsDataMap);

        // restrict to the loaded ref samples
        for(SampleTraitsData traits : traitsDataMap.values())
        {
            if(mSampleDataCache.hasRefSample(traits.SampleId))
            {
                sampleTraitsData.put(traits.SampleId, traits);
                assignSampleTraitsData(traits);
            }
        }
    }
}
