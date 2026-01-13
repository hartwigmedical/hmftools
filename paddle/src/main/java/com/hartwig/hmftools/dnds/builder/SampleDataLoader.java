package com.hartwig.hmftools.dnds.builder;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;
import static com.hartwig.hmftools.dnds.DndsCommon.MAX_REPEAT_COUNT;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.byteToBoolean;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import static org.jooq.impl.DSL.count;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.HotspotType;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.dnds.SampleMutationalLoad;
import com.hartwig.hmftools.dnds.SomaticVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jooq.Record;
import org.jooq.Record10;
import org.jooq.Record3;
import org.jooq.Result;

public class SampleDataLoader
{
    private final DatabaseAccess mDbAccess;
    private final String mPurpleDir;

    public SampleDataLoader(final String purpleDir, final DatabaseAccess dbAccess)
    {
        mDbAccess = dbAccess;
        mPurpleDir = purpleDir;
    }

    public List<SomaticVariant> loadVariants(final String sampleId)
    {
        if(mDbAccess != null)
        {
            return readDndsVariants(sampleId, MAX_REPEAT_COUNT);
        }
        else
        {
            DN_LOGGER.error("currently unsupported");
            return Collections.emptyList();
        }
    }

    public SampleMutationalLoad calcSampleMutationalLoad(final String sampleId)
    {
        if(mDbAccess != null)
        {
            return readSampleMutationalLoad(sampleId);
        }
        else
        {
            DN_LOGGER.error("currently unsupported");
            return null;
        }
    }

    private SampleMutationalLoad readSampleMutationalLoad(final String sample)
    {
        Result<Record3<Byte, String, Integer>> result = mDbAccess.context()
                .select(SOMATICVARIANT.BIALLELIC, SOMATICVARIANT.TYPE, count())
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.SAMPLEID.eq(sample))
                .and(SOMATICVARIANT.TYPE.in(VariantType.INDEL.toString(), VariantType.SNP.toString()))
                .and(SOMATICVARIANT.FILTER.eq("PASS"))
                .groupBy(SOMATICVARIANT.BIALLELIC, SOMATICVARIANT.TYPE)
                .fetch();

        int snvBiallelic = 0;
        int snvNonBiallelic = 0;
        int indelBiallelic = 0;
        int indelNonBiallelic = 0;

        for(Record3<Byte, String, Integer> record : result)
        {
            boolean isBiallelic = byteToBoolean(record.value1());
            VariantType type = VariantType.valueOf(record.value2());
            int count = record.value3();
            if(isBiallelic && type == VariantType.INDEL)
            {
                indelBiallelic = count;
            }
            if(isBiallelic && type == VariantType.SNP)
            {
                snvBiallelic = count;
            }
            if(!isBiallelic && type == VariantType.INDEL)
            {
                indelNonBiallelic = count;
            }
            if(!isBiallelic && type == VariantType.SNP)
            {
                snvNonBiallelic = count;
            }
        }

        PurityContext purityContext = mDbAccess.readPurityContext(sample);

        double purity = 0;

        if(purityContext == null)
        {
            DN_LOGGER.warn("sample({}) missing purity info", sample);
        }
        else
        {
            purity = purityContext.bestFit().purity();
        }

        return new SampleMutationalLoad(purity, indelBiallelic, indelNonBiallelic, snvBiallelic, snvNonBiallelic);
    }

    private List<SomaticVariant> readDndsVariants(final String sample, int maxRepeatCount)
    {
        List<SomaticVariant> variants = Lists.newArrayList();

        Result<Record10<String,Integer,String,String,String,Byte,String,String,String,Integer>> result = mDbAccess.context()
                .select(SOMATICVARIANT.CHROMOSOME, SOMATICVARIANT.POSITION, SOMATICVARIANT.REF, SOMATICVARIANT.ALT, SOMATICVARIANT.GENE,
                        SOMATICVARIANT.BIALLELIC, SOMATICVARIANT.HOTSPOT, SOMATICVARIANT.WORSTCODINGEFFECT,
                        SOMATICVARIANT.CANONICALCODINGEFFECT, SOMATICVARIANT.REPEATCOUNT)
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.SAMPLEID.eq(sample))
                .and(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.GENE.ne(""))
                .and(SOMATICVARIANT.REPEATCOUNT.lessOrEqual(maxRepeatCount))
                .and(SOMATICVARIANT.TYPE.in(VariantType.INDEL.toString(), VariantType.SNP.toString()))
                .fetch();

        for(Record record : result)
        {
            SomaticVariant variant = new SomaticVariant(
                    record.getValue(SOMATICVARIANT.CHROMOSOME),
                    record.getValue(SOMATICVARIANT.POSITION),
                    record.getValue(SOMATICVARIANT.REF),
                    record.getValue(SOMATICVARIANT.ALT),
                    record.getValue(SOMATICVARIANT.GENE),
                    byteToBoolean(record.getValue(SOMATICVARIANT.BIALLELIC)),
                    HotspotType.valueOf(record.getValue(SOMATICVARIANT.HOTSPOT)).equals(HotspotType.HOTSPOT),
                    record.getValue(SOMATICVARIANT.WORSTCODINGEFFECT).isEmpty()
                            ? CodingEffect.UNDEFINED
                            : CodingEffect.valueOf(record.getValue(SOMATICVARIANT.WORSTCODINGEFFECT)),
                    record.getValue(SOMATICVARIANT.CANONICALCODINGEFFECT).isEmpty()
                            ? CodingEffect.UNDEFINED
                            : CodingEffect.valueOf(record.getValue(SOMATICVARIANT.CANONICALCODINGEFFECT)),
                    record.getValue(SOMATICVARIANT.REPEATCOUNT));

            variants.add(variant);
        }

        return variants;
    }


}
