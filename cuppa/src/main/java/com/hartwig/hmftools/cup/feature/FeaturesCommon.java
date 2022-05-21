package com.hartwig.hmftools.cup.feature;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.feature.FeatureType.AMP;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.AMP_CN;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_AMP;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.traits.SampleTraitsData;

import org.apache.commons.cli.Options;

public final class FeaturesCommon
{
    public static final String COMBINE_DRIVER_AMP = "combine_driver_amp";
    public static final String RESTRICT_DRIVER_AMP_GENES = "restrict_driver_amp_genes";
    public static final String MIN_AMP_MULTIPLE = "min_amp_multiple";

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(COMBINE_DRIVER_AMP, false, "Combine driver AMPs with other driver events");
        options.addOption(RESTRICT_DRIVER_AMP_GENES, false, "Restrict driver AMPs to specifc list of genes");
        options.addOption(MIN_AMP_MULTIPLE, true, "AMP CN multiple of ploidy for use as a feature");
    }

    public static final List<String> DRIVER_AMP_GENES = Lists.newArrayList(
            "MYC","CCND1","ERBB2","ZNF703","FGFR1","MDM2","AR","CCNE1","ZNF217","EGFR","KRAS","CDK4","TERT","NCOA2","MET","MCL1",
            "ZMIZ1","FOXA1","CRYBG1","GNAS","CD44","CCND2","MDM4","CDX2","CDK6","TYMS","SOX4","MECOM","VEGFA","IRS2","ADAM30","PIK3CA",
            "GATA6","CD274","KIT","MYCL","KLF5","ARID5B","RAF1");

    public static void convertAndFilterDriverAmps(final Map<String,List<SampleFeatureData>> sampleFeaturesMap, boolean restrictGenes)
    {
        for(Map.Entry<String,List<SampleFeatureData>> sampleEntry : sampleFeaturesMap.entrySet())
        {
            List<SampleFeatureData> features = sampleEntry.getValue();

            for(int i = 0; i < features.size(); ++i)
            {
                SampleFeatureData featureData = features.get(i);

                if(featureData.Type != DRIVER || !featureData.ExtraInfo.containsValue(DRIVER_TYPE_AMP))
                    continue;

                if(restrictGenes && !DRIVER_AMP_GENES.contains(featureData.Name))
                    continue;

                SampleFeatureData ampFeature = new SampleFeatureData(featureData.SampleId, featureData.Name, AMP, featureData.Likelihood);
                ampFeature.ExtraInfo.putAll(featureData.ExtraInfo);

                features.set(i, ampFeature);
            }
        }
    }

    public static void filterDriverAmps(
            final Map<String,List<SampleFeatureData>> sampleFeaturesMap, final Map<String, SampleTraitsData> sampleTraitsData,
            double minCnMultiple)
    {
        for(Map.Entry<String,List<SampleFeatureData>> sampleEntry : sampleFeaturesMap.entrySet())
        {
            String sampleId = sampleEntry.getKey();
            SampleTraitsData traitsData = sampleTraitsData.get(sampleId);

            if(traitsData == null)
            {
                CUP_LOGGER.error("sample({}) missing traits data", sampleId);
                return;
            }

            double samplePloidy = traitsData.Ploidy;

            List<SampleFeatureData> features = sampleEntry.getValue();

            int index = 0;
            while(index < features.size())
            {
                SampleFeatureData featureData = features.get(index);

                if(featureData.Type == AMP && featureData.ExtraInfo.containsKey(AMP_CN))
                {
                    double adjustedCopyNumber = Double.parseDouble(featureData.ExtraInfo.get(AMP_CN)) / samplePloidy;
                    if(adjustedCopyNumber < minCnMultiple)
                    {
                        features.remove(index);
                        continue;
                    }
                }

                ++index;
            }
        }
    }
}
