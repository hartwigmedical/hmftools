package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

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
    private static final String COL_COHORT_NAME = "CohortName";

    public int sampleCountInCohort(final List<String> sampleIds, final String cohortName)
    {
        return (int)sampleIds.stream().map(x -> SampleCohort.get(x)).filter(y -> y.equals(cohortName)).count();
    }

    public boolean sampleInCohort(final String sampleId, final String cohortName)
    {
        return cohortName.equals(SampleCohort.get(sampleId));
    }

    public void removeMissingSamples(final List<String> sampleIds)
    {
        for(String sampleId : sampleIds)
        {
            SampleIds.remove(sampleId);
            SampleCancerType.remove(sampleId);
            SampleCohort.remove(sampleId);

            for(List<String> cancerSamples : CancerTypeSamples.values())
            {
                cancerSamples.remove(sampleId);
            }
        }
    }

    private boolean loadData(final String inputFile)
    {
        try
        {
            final List<String> items = Files.readAllLines(new File(inputFile).toPath());

            if(items.isEmpty())
                return false;

            final Map<String,Integer> fieldIndexMap = createFieldsIndexMap(items.get(0), DELIMITER);
            items.remove(0);

            Integer cohortNameIndex = fieldIndexMap.get(COL_COHORT_NAME);

            for(String item : items)
            {
                final String[] data = item.split(DELIMITER);
                if(data.length < 2)
                    return false;

                final String sampleId = data[COL_SAMPLE_ID];
                final String cancerType = data[COL_CANCER_TYPE];

                SampleIds.add(sampleId);
                SampleCancerType.put(sampleId, cancerType);

                if(cohortNameIndex != null)
                {
                    final String cohortName = data[cohortNameIndex];

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
