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
    public final Map<String, List<String>> CancerTypeSamples;
    public final Map<String, String> SampleCancerType;
    public final Map<String, String> SampleMutationType;
    public final List<String> SampleIds;
    public final List<String> MutationTypes;
    private boolean mIsValid;

    public SampleDataCache(final String inputFile)
    {
        SampleIds = Lists.newArrayList();
        CancerTypeSamples = Maps.newHashMap();
        SampleCancerType = Maps.newHashMap();
        SampleMutationType = Maps.newHashMap();
        MutationTypes = Lists.newArrayList();

        mIsValid = loadData(inputFile);
    }

    public boolean isValid() { return mIsValid; }

    private static final int COL_SAMPLE_ID = 0;
    private static final int COL_CANCER_TYPE = 1;
    private static final int COL_MUTATION_TYPE = 2;
    private static final String NO_MUTATION_TYPE = "NONE";

    public int sampleCountWithMutation(final List<String> sampleIds)
    {
        return (int)sampleIds.stream().filter(x -> SampleMutationType.containsKey(x)).count();
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

                if(item.length() > COL_MUTATION_TYPE)
                {
                    final String mutationType = data[COL_MUTATION_TYPE];

                    if(!mutationType.equals(NO_MUTATION_TYPE))
                    {
                        SampleMutationType.put(sampleId, mutationType);

                        if(!MutationTypes.contains(mutationType))
                            MutationTypes.add(mutationType);
                    }
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
