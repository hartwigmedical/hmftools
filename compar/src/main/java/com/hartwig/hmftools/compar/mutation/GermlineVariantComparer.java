package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.mutation.VariantData.addComparerFields;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldCheck;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineVariantComparer extends ItemComparer
{
    public GermlineVariantComparer(final ComparConfig config, final Map<String,FieldCheck> fieldCheckMap)
    {
        super(config);

        addComparerFields(mFields, fieldCheckMap);
    }

    @Override
    public CategoryType category() { return GERMLINE_VARIANT; }

    @Override
    public List<String> displayFieldNames()
    {
        return VariantData.sharedFieldNames();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        String vcfFile = PurpleCommon.purpleGermlineVcfFile(fileSources.Purple, sampleId);

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        if(!vcfFileReader.fileValid())
        {
            CMP_LOGGER.warn("failed to read germline VCF file({})", vcfFile);
            return null;
        }

        try
        {
            for(VariantContext variantContext : vcfFileReader.iterator())
            {
                GermlineVariantData variant = GermlineVariantData.fromContext(
                        variantContext, sampleId, false, fileSources.Source, mConfig, mFields);

                if(mConfig.RestrictToDrivers && !mConfig.DriverGenes.contains(variant.Gene))
                    continue;

                comparableItems.add(variant);
            }

            CMP_LOGGER.debug("sample({}) loaded {} {} germline variants", sampleId, fileSources.Source, comparableItems.size());
        }
        catch(Exception e)
        {
            CMP_LOGGER.warn("failed to read germline VCF file({}): {}", vcfFile, e.toString());
            return null;
        }

        return comparableItems;
    }
}
