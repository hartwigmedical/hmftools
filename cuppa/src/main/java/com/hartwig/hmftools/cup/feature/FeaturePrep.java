package com.hartwig.hmftools.cup.feature;

import static com.hartwig.hmftools.common.cuppa.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromFile;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

public class FeaturePrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    public FeaturePrep(final PrepConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType categoryType() { return FEATURE; }

    private static final String AMP_TYPE = ".amp";
    private static final String MUTATION_TYPE = ".mut";
    private static final String INDEL_TYPE = ".indel";

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        LinkedHashMap<DataItem.Index, DataItem> dataItemsMap = new LinkedHashMap<>();

        final String linxDataDir = formSamplePath(mConfig.LinxDir, sampleId);
        final String purpleDataDir = formSamplePath(mConfig.PurpleDir, sampleId);
        final String virusDataDir = formSamplePath(mConfig.VirusDir, sampleId);

        try
        {
            // TODO: remove use of sample map once Cuppa becomes just a data-loader, or refactor
            Map<String,List<SampleFeatureData>> sampleFeaturesMap = new HashMap<>();

            if(!loadFeaturesFromFile(sampleId, linxDataDir, purpleDataDir, virusDataDir, sampleFeaturesMap))
                return null;

            /* Examples:

            Name,Type,Likelihood,ExtraInfo
            TP53,DRIVER,1.0,CHR=17
            PABPC1,DRIVER,1.0,CN=14.8;CHR=8;TYPE=AMP
            RAD51B,DRIVER,0.0,CHR=14
            HPV,VIRUS,1.0,
            TMPRSS2_ERG,FUSION,1.0,
            INDEL_SLC34A2,INDEL,1.0,
            */

            for(SampleFeatureData featureData : sampleFeaturesMap.get(sampleId))
            {
                String likelihood = String.valueOf(featureData.Likelihood);

                DataItem dataItem;
                switch(featureData.Type)
                {
                    case DRIVER:
                        dataItem = new DataItem(DNA, ItemType.DRIVER, featureData.Name + MUTATION_TYPE, likelihood);
                        break;

                    case AMP:
                        dataItem = new DataItem(DNA, ItemType.DRIVER, featureData.Name + AMP_TYPE, likelihood);
                        break;

                    case INDEL:
                        dataItem = new DataItem(DNA, ItemType.DRIVER, featureData.Name.replace("INDEL_", "") + INDEL_TYPE, likelihood);
                        break;

                    case FUSION:
                        dataItem = new DataItem(DNA, ItemType.FUSION, featureData.Name, likelihood);
                        break;

                    case VIRUS:
                        dataItem = new DataItem(DNA, ItemType.VIRUS, featureData.Name, likelihood);
                        break;

                    default:
                        throw new IllegalStateException("Invalid FeatureType");
                }

                // De-duplicate features by max value
                if(dataItemsMap.containsKey(dataItem.Index))
                {
                    float newDataItemValue = Float.parseFloat(dataItem.Value);
                    float existingDataItemValue = Float.parseFloat(dataItemsMap.get(dataItem.Index).Value);

                    if(existingDataItemValue >= newDataItemValue)
                    {
                        // Keep old value if it is higher and discard new value
                        continue;
                    }
                }

                dataItemsMap.put(dataItem.Index, dataItem);
            }

            return new ArrayList<>(dataItemsMap.values());
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) - failed to extract driver, fusion, and virus features: {}", sampleId, e.toString());
            return null;
        }
    }
}
