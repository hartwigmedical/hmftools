package com.hartwig.hmftools.cup.feature;

import static com.hartwig.hmftools.cup.feature.FeatureType.AMP;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_AMP;

import java.util.List;
import java.util.Map;

public final class FeaturesCommon
{
    public static void convertDriverAmps(final Map<String,List<SampleFeatureData>> sampleFeaturesMap)
    {
        for(Map.Entry<String,List<SampleFeatureData>> sampleEntry : sampleFeaturesMap.entrySet())
        {
            List<SampleFeatureData> features = sampleEntry.getValue();

            for(int i = 0; i < features.size(); ++i)
            {
                SampleFeatureData featureData = features.get(i);

                if(featureData.Type != DRIVER || !featureData.ExtraInfo.containsValue(DRIVER_TYPE_AMP))
                    continue;

                SampleFeatureData ampFeature = new SampleFeatureData(featureData.SampleId, featureData.Name, AMP, featureData.Likelihood);
                ampFeature.ExtraInfo.putAll(featureData.ExtraInfo);

                features.set(i, ampFeature);
            }
        }
    }
}
