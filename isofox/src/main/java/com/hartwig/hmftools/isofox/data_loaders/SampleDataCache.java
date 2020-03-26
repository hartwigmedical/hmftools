package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class SampleDataCache
{
    public final Map<String, List<String>> CancerTypeSamples; // cancer type to list of samples
    public final Map<String, String> SampleCancerType; // sample to cancer type
    public final Map<String, String> SampleCohort; // sample to cohort name
    public final List<String> SampleIds;
    public final List<String> CohortNames;
    private boolean mIsValid;

    public final static String COHORT_A = "CohortA"; // 'cohort A' being evaluated, may be arbitrary
    public final static String COHORT_B = "CohortB"; // 'cohort A' being evaluated, may be arbitrary

    public SampleDataCache(final String inputFile)
    {
        SampleIds = Lists.newArrayList();
        CancerTypeSamples = Maps.newHashMap();
        SampleCancerType = Maps.newHashMap();
        SampleCohort = Maps.newHashMap();
        CohortNames = Lists.newArrayList();

        mIsValid = loadData(inputFile);
    }

    public boolean isValid() { return mIsValid; }

    private static final int COL_SAMPLE_ID = 0;
    private static final int COL_CANCER_TYPE = 1;
    private static final int COL_COHORT_NAME = 2;

    public int sampleCountInCohort(final List<String> sampleIds, final String cohortName)
    {
        return (int)sampleIds.stream().map(x -> SampleCohort.get(x)).filter(y -> y.equals(cohortName)).count();
    }

    public boolean sampleInCohort(final String sampleId, final String cohortName)
    {
        return cohortName.equals(SampleCohort.get(sampleId));
    }

    private boolean loadData(final String inputFile)
    {
        try
        {
            final List<String> items = Files.readAllLines(new File(inputFile).toPath());

            for(String item : items)
            {
                if(item.contains("SampleId"))
                    continue;

                final String[] data = item.split(",");
                if(data.length < 2)
                    return false;

                final String sampleId = data[COL_SAMPLE_ID];
                final String cancerType = data[COL_CANCER_TYPE];

                SampleIds.add(sampleId);
                SampleCancerType.put(sampleId, cancerType);

                if(item.length() > COL_COHORT_NAME)
                {
                    final String cohortName = data[COL_COHORT_NAME];

                    SampleCohort.put(sampleId, cohortName);

                    if(!CohortNames.contains(cohortName))
                        CohortNames.add(cohortName);
                }

                List<String> sampleIds = CancerTypeSamples.get(cancerType);
                if(sampleIds == null)
                {
                    sampleIds = Lists.newArrayList();
                    CancerTypeSamples.put(cancerType, sampleIds);
                }

                sampleIds.add(sampleId);
            }

            ISF_LOGGER.info("load {} samples, {} cancer types from file({})", SampleIds.size(), CancerTypeSamples.size(), inputFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load cancer-type sample file({}): {}", inputFile, e.toString());
            return false;
        }

        return true;
    }

}
