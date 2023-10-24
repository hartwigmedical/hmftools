package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.cuppa.CategoryType.SV;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.purpleSvFile;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.extractSvCounts;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.cuppa.SvDataType;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

public class StructuralVariantPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    public StructuralVariantPrep(final PrepConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType categoryType() { return SV; }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        List<DataItem> dataItems = Lists.newArrayList();

        final String purpleDataDir = mConfig.getPurpleDataDir(sampleId);

        try
        {
            final String svVcfFile = purpleSvFile(mConfig.PurpleDir, sampleId);
            final String linxDataDir = mConfig.getLinxDataDir(sampleId);
            final String clusterFile = LinxCluster.generateFilename(linxDataDir, sampleId, false);

            final List<StructuralVariantData> svDataList = Lists.newArrayList();
            final List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(svVcfFile, new AlwaysPassFilter());
            final List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            int svId = 0;
            for (EnrichedStructuralVariant var : enrichedVariants)
            {
                svDataList.add(convertSvData(var, svId++));
            }

            final List<LinxCluster> clusterList = LinxCluster.read(clusterFile);

            final int[] svDataCounts = extractSvCounts(svDataList, clusterList);

            for(SvDataType type : SvDataType.values())
            {
                dataItems.add(new DataItem(DNA, ItemType.SV_COUNT, type.toString(), String.valueOf(svDataCounts[type.ordinal()])));
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
