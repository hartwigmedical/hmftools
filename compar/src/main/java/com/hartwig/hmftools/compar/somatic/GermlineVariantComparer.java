package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_GERMLINE_VCF_SUFFIX;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.compar.Category.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.DiffFunctions.diffsStr;
import static com.hartwig.hmftools.compar.Mismatch.commonCsv;
import static com.hartwig.hmftools.compar.somatic.VariantCommon.FLD_QUAL;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GERMLINEVARIANT;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.SomaticVariantDAO;

import org.jooq.Record;
import org.jooq.Result;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.vcf.VCFCodec;

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
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public boolean hasDetailedOutput() { return true; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        // same as somatic
        thresholds.addFieldThreshold(FLD_QUAL, 20, 0.2);
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
            SomaticVariant variant = SomaticVariantDAO.buildFromRecord(record);
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

        SomaticVariantFactory variantFactory = new SomaticVariantFactory(filter);

        String vcfFile = !fileSources.GermlineVcf.isEmpty() ?
                fileSources.GermlineVcf : fileSources.Purple + sampleId + PURPLE_GERMLINE_VCF_SUFFIX;

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);

            for(VariantContext variantContext : reader.iterator())
            {
                if(filter.test(variantContext))
                {
                    final SomaticVariant variant = variantFactory.createVariant(sampleId, variantContext).orElse(null);

                    if(variant == null)
                        continue;

                    comparableItems.add(new GermlineVariantData(variant));
                }
            }

            CMP_LOGGER.debug("sample({}) loaded {} germline variants", sampleId, comparableItems.size());
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to read germline VCF file({}): {}", vcfFile, e.toString());
        }

        return comparableItems;
    }

    @Override
    public String outputHeader()
    {
        return "Category,MismatchType,Key,RefQual,NewQual,RefTier,NewTier,RefTotalReads,NewTotalReads,RefAlleleReads,NewAlleleReads,Differences";
    }

    @Override
    public String mismatchOutput(final Mismatch mismatch)
    {
        final GermlineVariantData refVar = mismatch.RefItem != null ? (GermlineVariantData)mismatch.RefItem : null;
        final GermlineVariantData otherVar = mismatch.NewItem != null ? (GermlineVariantData)mismatch.NewItem : null;

        return String.format("%s,%.0f,%.0f,%s,%s,%d,%d,%d,%d,%s",
                commonCsv(mismatch),
                refVar != null ? refVar.Variant.qual() : 0, otherVar != null ? otherVar.Variant.qual() : 0,
                refVar != null ? refVar.Variant.tier() : "", otherVar != null ? otherVar.Variant.tier() : "",
                refVar != null ? refVar.Variant.totalReadCount() : 0, otherVar != null ? otherVar.Variant.totalReadCount() : 0,
                refVar != null ? refVar.Variant.alleleReadCount() : 0, otherVar != null ? otherVar.Variant.alleleReadCount() : 0,
                diffsStr(mismatch.DiffValues));
    }

}
