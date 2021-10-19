package com.hartwig.hmftools.neo.epitope;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;

import java.util.List;

import org.apache.commons.compress.utils.Lists;

public class SampleData
{
    public final String Id;
    public final String CancerType;
    public final List<String> HlaTypes;

    public SampleData(final String id, final String cancerType, final List<String> hlaTypes)
    {
        Id = id;
        CancerType = cancerType;
        HlaTypes = hlaTypes;
    }

    public static SampleData fromCsv(final String data)
    {
        final String[] items = data.split(DELIMITER, -1);

        if(items.length < 2)
            return null;

        String sampleId = items[0];
        String cancerType = items[1];

        final List<String> hlaTypes = Lists.newArrayList();

        if(items.length == 8)
        {
            for(int i = 2; i < 8; ++i)
            {
                hlaTypes.add(items[i]);
            }
        }

        return new SampleData(sampleId, cancerType, hlaTypes);
    }
}
