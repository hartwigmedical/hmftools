package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.compar.driver.DriverData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.SomaticVariantDAO;

import org.apache.commons.compress.utils.Lists;
import org.jooq.Record;
import org.jooq.Result;

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

        final List<List<ComparableItem>> sourceDrivers = Lists.newArrayList();

        for(String sourceName : mConfig.DbSourceNames)
        {
            sourceDrivers.add(getSampleVariants(sampleId, mConfig.DbConnections.get(sourceName)));
        }

        for(int i = 0; i < mConfig.DbSourceNames.size() - 1; ++i)
        {
            final String source1 = mConfig.DbSourceNames.get(i);

            for(int j = i + 1; j < mConfig.DbSourceNames.size(); ++j)
            {
                final String source2 = mConfig.DbSourceNames.get(j);

                CommonUtils.compareItems(sampleId, mismatches, matchLevel, source1, source2, sourceDrivers.get(i), sourceDrivers.get(j));
            }
        }
    }

    private List<ComparableItem> getSampleVariants(final String sampleId, final DatabaseAccess dbAccess)
    {
        Result<Record> result = dbAccess.context().select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.GENE.isNotNull())
                .and(SOMATICVARIANT.WORSTCODINGEFFECT.in(NONSENSE_OR_FRAMESHIFT.toString(), SPLICE.toString(), MISSENSE.toString()))
                .fetch();

        final List<ComparableItem> variants = Lists.newArrayList();
        for (Record record : result)
        {
            SomaticVariant variant = SomaticVariantDAO.buildFromRecord(record);
            variants.add(new SomaticVariantData(variant));
        }

        return variants;
    }


}
