package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.compar.driver.DriverData;
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

    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        final MatchLevel matchLevel = mConfig.Categories.get(SOMATIC_VARIANT);

        final List<List<ComparableItem>> sourceItems = Lists.newArrayList();

        for(String sourceName : mConfig.SourceNames)
        {
            String sourceSampleId = mConfig.sourceSampleId(sourceName, sampleId);

            if(!mConfig.DbConnections.isEmpty())
                sourceItems.add(loadFromDb(sourceSampleId, mConfig.DbConnections.get(sourceName)));
            else
                sourceItems.add(loadFromFile(sourceSampleId, mConfig.FileSources.get(sourceName)));
        }

        for(int i = 0; i < mConfig.SourceNames.size() - 1; ++i)
        {
            final String source1 = mConfig.SourceNames.get(i);

            for(int j = i + 1; j < mConfig.SourceNames.size(); ++j)
            {
                final String source2 = mConfig.SourceNames.get(j);

                CommonUtils.compareItems(sampleId, mismatches, matchLevel, source1, source2, sourceItems.get(i), sourceItems.get(j));
            }
        }
    }

    private List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        Result<Record> result = dbAccess.context().select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.GENE.isNotNull())
                .and(SOMATICVARIANT.CANONICALCODINGEFFECT.in(NONSENSE_OR_FRAMESHIFT.toString(), SPLICE.toString(), MISSENSE.toString()))
                .fetch();

        final List<ComparableItem> variants = Lists.newArrayList();
        for (Record record : result)
        {
            SomaticVariant variant = SomaticVariantDAO.buildFromRecord(record);
            variants.add(new SomaticVariantData(variant));
        }

        return variants;
    }

    private List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        SomaticVariantFactory variantFactory = new SomaticVariantFactory(filter);

        String vcfFile = fileSources.Purple + sampleId + PURPLE_SOMATIC_VCF_SUFFIX;

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);

            for(VariantContext variant : reader.iterator())
            {
                if(filter.test(variant))
                {
                    final SomaticVariant somaticVariant = variantFactory.createVariant(sampleId, variant).orElse(null);

                    if(somaticVariant == null)
                        continue;

                    comparableItems.add(new SomaticVariantData(somaticVariant));
                }
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to read Purple somatic VCF file({}): {}", vcfFile, e.toString());
        }

        return comparableItems;
    }
}
