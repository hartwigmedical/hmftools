package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
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

    /*
    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new BooleanField(FLD_REPORTED, i -> ((GermlineVariantData) i).Variant.reported(), true),
                new StringField(FLD_HOTSPOT, i -> ((GermlineVariantData) i).Variant.hotspot().toString(), true),
                new StringField(FLD_TIER, i -> ((GermlineVariantData) i).Variant.tier().toString(), true),
                new StringField(FLD_GENE, i -> ((GermlineVariantData) i).Variant.gene(), true),
                new StringField(FLD_CANON_EFFECT, i -> ((GermlineVariantData) i).Variant.canonicalEffect(), true),
                new StringField(FLD_CODING_EFFECT, i -> ((GermlineVariantData) i).Variant.canonicalCodingEffect().toString(), true),
                new StringField(FLD_HGVS_CODING, i -> ((GermlineVariantData) i).Variant.canonicalHgvsCodingImpact(), true),
                new StringField(FLD_HGVS_PROTEIN, i -> ((GermlineVariantData) i).Variant.canonicalHgvsProteinImpact(), true),
                new StringField(FLD_OTHER_REPORTED, i -> ((GermlineVariantData) i).Variant.otherReportedEffects(), true),
                new IntField(FLD_QUAL, i -> (int) ((GermlineVariantData) i).Variant.qual(), true,
                        50.0, 0.2),
                new DoubleField(FLD_VARIANT_COPY_NUMBER, i -> ((GermlineVariantData) i).Variant.variantCopyNumber(),
                        true, 0.3, 0.3, "%.2f"),
                new DoubleField(FLD_PURITY_ADJUSTED_VAF, i -> ((GermlineVariantData) i).Variant.adjustedVAF(),
                        true, 0.2, null, "%.2f"),
                new IntField(FLD_TUMOR_SUPPORTING_READ_COUNT,
                        i -> ((GermlineVariantData) i).Variant.allelicDepth().AlleleReadCount, true,
                        1.0, 0.2),
                new IntField(FLD_TUMOR_TOTAL_READ_COUNT,
                        i -> ((GermlineVariantData) i).Variant.allelicDepth().TotalReadCount, true,
                        1.0, 0.2),
                new StringListField(FLD_FILTER, i -> ((GermlineVariantData) i).Filters.stream().sorted().toList(), true)
        );
    }
    */

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
                if(variantContext.isFiltered())
                    continue;

                GermlineVariantData variant = GermlineVariantData.fromContext(
                        variantContext, sampleId, false, fileSources.Source, mConfig, mFields);

                if(mConfig.RestrictToDrivers && !mConfig.DriverGenes.contains(variant.Gene))
                    continue;

                comparableItems.add(variant);
            }

            /*
            List<SmallVariant> variants = SmallVariantFactory.loadVariants(sampleId, vcfFile);
            for(SmallVariant variant : variants)
            {
                if(!variant.isFiltered())
                {
                    BasePosition comparisonPosition = determineComparisonGenomePosition(
                            variant.chromosome(), variant.position(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                    GermlineVariantData variant = GermlineVariantData.fromContext(
                            variantContext, sampleId, false, fileSources.Source, mConfig, mFields);

                    comparableItems.add(new GermlineVariantData(variant, comparisonPosition, mFields));
                }
            }
            */

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
