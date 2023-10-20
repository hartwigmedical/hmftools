package com.hartwig.hmftools.cup.feature;

import static com.hartwig.hmftools.common.cuppa.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;
import static com.hartwig.hmftools.cup.common.CupCalcs.BOOL_STR_TRUE;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromFile;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
    private static final String MUTATION_TYPE = ".mutation";
    private static final String INDEL_TYPE = ".indel";

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        List<DataItem> dataItems = Lists.newArrayList();

        final String linxDataDir = formSamplePath(mConfig.LinxDir, sampleId);
        final String purpleDataDir = formSamplePath(mConfig.PurpleDir, sampleId);
        final String virusDataDir = formSamplePath(mConfig.VirusDir, sampleId);

        try
        {
            // TODO: remove use of sample map once Cuppa becomes just a data-loader, or refactor
            Map<String,List<SampleFeatureData>> sampleFeaturesMap = Maps.newHashMap();

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
                switch(featureData.Type)
                {
                    case DRIVER:
                    {
                        dataItems.add(new DataItem(DNA, ItemType.DRIVER, featureData.Name + MUTATION_TYPE, BOOL_STR_TRUE));
                        break;
                    }

                    case AMP:
                    {
                        dataItems.add(new DataItem(
                                DNA, ItemType.DRIVER, featureData.Name + AMP_TYPE, String.valueOf(featureData.Likelihood)));
                        break;
                    }

                    case INDEL:
                    {
                        dataItems.add(new DataItem(
                                DNA, ItemType.DRIVER, featureData.Name + INDEL_TYPE, String.valueOf(featureData.Likelihood)));
                        break;
                    }

                    case FUSION:
                    {
                        dataItems.add(new DataItem(DNA, ItemType.FUSION, featureData.Name, BOOL_STR_TRUE));
                        break;
                    }

                    case VIRUS:
                    {
                        dataItems.add(new DataItem(DNA, ItemType.VIRUS, featureData.Name, BOOL_STR_TRUE));
                        break;
                    }
                }
            }

            return dataItems;
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) sample traits - failed to load purity file from dir{}): {}",
                    sampleId, purpleDataDir, e.toString());

            return null;
        }
    }
}
