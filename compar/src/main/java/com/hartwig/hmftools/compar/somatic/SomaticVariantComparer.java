package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.DiffFunctions.diffsStr;
import static com.hartwig.hmftools.compar.Mismatch.commonCsv;
import static com.hartwig.hmftools.compar.Mismatch.commonHeader;
import static com.hartwig.hmftools.compar.somatic.SomaticVariantData.FLD_QUAL;
import static com.hartwig.hmftools.compar.somatic.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;

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

public class SomaticVariantComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public SomaticVariantComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return SOMATIC_VARIANT; }

    @Override
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_QUAL, 20, 0.2);
        thresholds.addFieldThreshold(FLD_SUBCLONAL_LIKELIHOOD, 0.6, 0);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        Result<Record> result = dbAccess.context().select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .fetch();

        final List<ComparableItem> variants = Lists.newArrayList();
        for (Record record : result)
        {
            SomaticVariant variant = SomaticVariantDAO.buildFromRecord(record);
            variants.add(new SomaticVariantData(variant));
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

        // use the Purple suffix if not specified
        String vcfFile = !fileSources.SomaticVcf.isEmpty() ?
                fileSources.SomaticVcf : fileSources.Purple + sampleId + PURPLE_SOMATIC_VCF_SUFFIX;

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);

            for(VariantContext variant : reader.iterator())
            {
                if(filter.test(variant))
                {
                    final SomaticVariant somaticVariant = variantFactory.createVariant(sampleId, variant).orElse(null);
                    comparableItems.add(new SomaticVariantData(somaticVariant));
                }
            }

            CMP_LOGGER.debug("sample({}) loaded {} somatic variants", sampleId, comparableItems.size());
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to read somatic VCF file({}): {}", vcfFile, e.toString());
        }

        return comparableItems;
    }

    @Override
    public String outputHeader()
    {
        return commonHeader() + ",RefFilter,RefQual,OtherQual,RefTier,OtherTier,RefTotalReads,NewTotalReads"
                + ",RefAlleleReads,NewAlleleReads,RefLPS,NewLPS,Differences";
    }

    @Override
    public String mismatchOutput(final Mismatch mismatch)
    {
        final SomaticVariantData refVar = mismatch.RefItem != null ? (SomaticVariantData)mismatch.RefItem : null;
        final SomaticVariantData newVar = mismatch.NewItem != null ? (SomaticVariantData)mismatch.NewItem : null;

        return String.format("%s,%s,%.0f,%.0f,%s,%s,%d,%d,%d,%d,%s,%s,%s",
                commonCsv(mismatch), refVar != null ? refVar.Variant.filter() : newVar.Variant.filter(),
                refVar != null ? refVar.Variant.qual() : 0, newVar != null ? newVar.Variant.qual() : 0,
                refVar != null ? refVar.Variant.tier() : "", newVar != null ? newVar.Variant.tier() : "",
                refVar != null ? refVar.Variant.totalReadCount() : 0, newVar != null ? newVar.Variant.totalReadCount() : 0,
                refVar != null ? refVar.Variant.alleleReadCount() : 0, newVar != null ? newVar.Variant.alleleReadCount() : 0,
                refVar != null ? refVar.Variant.localPhaseSetsStr() : "", newVar != null ? newVar.Variant.localPhaseSetsStr() : "",
                diffsStr(mismatch.DiffValues));
    }

}
