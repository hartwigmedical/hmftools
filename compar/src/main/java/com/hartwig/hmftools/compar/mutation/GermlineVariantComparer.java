package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
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
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
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
    public void registerThresholds(final FieldConfig fieldConfig)
    {
        VariantCommon.registerThresholds(category(), fieldConfig);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return VariantCommon.comparedFieldNames();
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
