package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.compar.Category.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GERMLINEVARIANT;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.GermlineVariantDAO;

import org.jooq.Record;
import org.jooq.Result;

import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class GermlineVariantComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public GermlineVariantComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return GERMLINE_VARIANT; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        // same as somatic
        thresholds.addFieldThreshold(FLD_QUAL, 20, 0.2);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        List<String> fieldNames = VariantCommon.comparedFieldNames();
        return fieldNames;
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        Result<Record> result = dbAccess.context().select()
                .from(GERMLINEVARIANT)
                .where(GERMLINEVARIANT.FILTER.eq(PASS_FILTER))
                .and(GERMLINEVARIANT.SAMPLEID.eq(sampleId))
                .fetch();

        final List<ComparableItem> variants = Lists.newArrayList();
        for (Record record : result)
        {
            GermlineVariant variant = GermlineVariantDAO.buildFromRecord(record);
            variants.add(new GermlineVariantData(variant));
        }

        return variants;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        String vcfFile = PurpleCommon.purpleGermlineVcfFile(fileSources.Purple, sampleId);

        try
        {
            List<GermlineVariant> variants = GermlineVariantFactory.fromVCFFile(sampleId, vcfFile);
            variants.forEach(x -> comparableItems.add(new GermlineVariantData(x)));

            CMP_LOGGER.debug("sample({}) loaded {} germline variants", sampleId, comparableItems.size());
        }
        catch(Exception e)
        {
            CMP_LOGGER.warn("failed to read germline VCF file({}): {}", vcfFile, e.toString());
            return null;
        }

        return comparableItems;
    }
}
