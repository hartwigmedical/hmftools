package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_FILTER;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_CANON_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_CODING_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_GENE;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HGVS_CODING;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HGVS_PROTEIN;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HOTSPOT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_OTHER_REPORTED;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_PURITY_ADJUSTED_VAF;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TIER;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TUMOR_SUPPORTING_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TUMOR_TOTAL_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_VARIANT_COPY_NUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GERMLINEVARIANT;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.variant.SmallVariantFactory;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.BooleanField;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.IntField;
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.compar.common.field.StringListField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.GermlineVariantDAO;

import org.jooq.Record;
import org.jooq.Result;

public class GermlineVariantComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public GermlineVariantComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return GERMLINE_VARIANT; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<Field> fields(MatchLevel matchLevel)
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
                new IntField(FLD_QUAL, i -> (int) ((GermlineVariantData) i).Variant.qual(), true, 50.0, 0.2),
                new DoubleField(FLD_VARIANT_COPY_NUMBER, i -> ((GermlineVariantData) i).Variant.variantCopyNumber(), true, 0.3, 0.3, "%.2f"),
                new DoubleField(FLD_PURITY_ADJUSTED_VAF, i -> ((GermlineVariantData) i).Variant.adjustedVAF(), true, 0.2, null, "%.2f"),
                new IntField(FLD_TUMOR_SUPPORTING_READ_COUNT,
                        i -> ((GermlineVariantData) i).Variant.allelicDepth().AlleleReadCount, true, 1.0, 0.2),
                new IntField(FLD_TUMOR_TOTAL_READ_COUNT,
                        i -> ((GermlineVariantData) i).Variant.allelicDepth().TotalReadCount, true, 1.0, 0.2),
                new StringListField(FLD_FILTER, i -> ((GermlineVariantData) i).Filters.stream().sorted().toList(), true)
        );
    }

    @Override
    public List<String> displayFieldNames()
    {
        return VariantCommon.sharedDisplayFieldNames();
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        Result<Record> result = dbAccess.context().select()
                .from(GERMLINEVARIANT)
                .where(GERMLINEVARIANT.FILTER.eq(PASS_FILTER))
                .and(GERMLINEVARIANT.SAMPLEID.eq(sampleId))
                .fetch();

        final List<ComparableItem> variants = Lists.newArrayList();
        for (Record record : result)
        {
            SmallVariant variant = GermlineVariantDAO.buildFromRecord(record);
            BasePosition comparisonPosition = determineComparisonGenomePosition(
                    variant.chromosome(), variant.position(), sourceType, mConfig.RequiresLiftover, mConfig.LiftoverCache);
            variants.add(new GermlineVariantData(variant, comparisonPosition));
        }

        return variants;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        String vcfFile = PurpleCommon.purpleGermlineVcfFile(fileSources.Purple, sampleId);

        try
        {
            List<SmallVariant> variants = SmallVariantFactory.loadVariants(sampleId, vcfFile);
            for(SmallVariant variant : variants)
            {
                if(!variant.isFiltered())
                {
                    BasePosition comparisonPosition = determineComparisonGenomePosition(
                            variant.chromosome(), variant.position(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);
                    comparableItems.add(new GermlineVariantData(variant, comparisonPosition));
                }
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
