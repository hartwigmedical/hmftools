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
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.HotspotType;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.dnds.SampleMutationalLoad;
import com.hartwig.hmftools.dnds.SomaticVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jooq.Record;
import org.jooq.Record11;
import org.jooq.Record3;
import org.jooq.Result;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.variantcontext.VariantContext;

public class SampleDataLoader
{
    private final DatabaseAccess mDbAccess;
    private final String mPurpleDir;

    public SampleDataLoader(final String purpleDir, final DatabaseAccess dbAccess)
    {
        mDbAccess = dbAccess;
        mPurpleDir = purpleDir;
    }

    public static class SampleData
    {
        public final List<SomaticVariant> Variants;
        public final SampleMutationalLoad MutationalLoad;

        public SampleData(final List<SomaticVariant> variants, final SampleMutationalLoad mutationalLoad)
        {
            Variants = variants;
            MutationalLoad = mutationalLoad;
        }
    }

    public SampleData loadSampleData(final String sampleId)
    {
        if(mDbAccess != null)
        {
            return new SampleData(readDndsVariants(sampleId, MAX_REPEAT_COUNT), readSampleMutationalLoad(sampleId));
        }
        else if(mPurpleDir != null)
        {
            return readSampleDataFromVcf(sampleId, MAX_REPEAT_COUNT);
        }
        else
        {
            DN_LOGGER.error("sample({}) neither database nor purple directory configured", sampleId);
            return new SampleData(Collections.emptyList(), null);
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

        return new SampleMutationalLoad(snvBiallelic, snvNonBiallelic, indelBiallelic, indelNonBiallelic);
    }

    private List<SomaticVariant> readDndsVariants(final String sample, int maxRepeatCount)
    {
        List<SomaticVariant> variants = Lists.newArrayList();

        Result<Record11<String,Integer,String,String,String,Byte,String,String,String,Byte,Integer>> result = mDbAccess.context()
                .select(SOMATICVARIANT.CHROMOSOME, SOMATICVARIANT.POSITION, SOMATICVARIANT.REF, SOMATICVARIANT.ALT, SOMATICVARIANT.GENE,
                        SOMATICVARIANT.BIALLELIC, SOMATICVARIANT.HOTSPOT, SOMATICVARIANT.WORSTCODINGEFFECT,
                        SOMATICVARIANT.CANONICALCODINGEFFECT, SOMATICVARIANT.SPLICEREGION, SOMATICVARIANT.REPEATCOUNT)
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
                    byteToBoolean(record.getValue(SOMATICVARIANT.SPLICEREGION)),
                    record.getValue(SOMATICVARIANT.REPEATCOUNT));

            variants.add(variant);
        }

        return variants;
    }

    private SampleData readSampleDataFromVcf(final String sample, final int maxRepeatCount)
    {
        String vcfFile = PurpleCommon.purpleSomaticVcfFile(mPurpleDir, sample);
        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        if(!vcfFileReader.fileValid())
        {
            DN_LOGGER.warn("sample({}) purple vcf file({}) not found, skipping", sample, vcfFile);
            return null;
        }

        int snvBiallelic = 0;
        int snvNonBiallelic = 0;
        int indelBiallelic = 0;
        int indelNonBiallelic = 0;

        List<SomaticVariant> variants = Lists.newArrayList();

        try(vcfFileReader; CloseableTribbleIterator<VariantContext> iterator = vcfFileReader.iterator())
        {
            for(VariantContext context : iterator)
            {
                VariantContextDecorator decorator = new VariantContextDecorator(context);

                if(!decorator.isPass())
                    continue;

                VariantType type = decorator.type();
                if(type != VariantType.SNP && type != VariantType.INDEL)
                    continue;

                boolean isBiallelic = decorator.biallelic();

                if(type == VariantType.SNP)
                {
                    if(isBiallelic)
                        ++snvBiallelic;
                    else
                        ++snvNonBiallelic;
                }
                else
                {
                    if(isBiallelic)
                        ++indelBiallelic;
                    else
                        ++indelNonBiallelic;
                }

                String gene = decorator.gene();
                if(gene.isEmpty())
                    continue;

                VariantImpact impact = decorator.variantImpact();
                CodingEffect codingEffect = impact.CanonicalCodingEffect;
                boolean hasNoCodingEffect = codingEffect == CodingEffect.NONE || codingEffect == CodingEffect.UNDEFINED;
                if(hasNoCodingEffect && !impact.CanonicalSpliceRegion)
                    continue;

                int repeatCount = decorator.repeatCount();
                if(repeatCount > maxRepeatCount)
                    continue;

                SomaticVariant variant = new SomaticVariant(
                        decorator.chromosome(),
                        decorator.position(),
                        decorator.ref(),
                        decorator.alt(),
                        gene,
                        isBiallelic,
                        decorator.isHotspot(),
                        impact.WorstCodingEffect,
                        codingEffect,
                        impact.CanonicalSpliceRegion,
                        repeatCount
                );

                variants.add(variant);
            }
        }

        SampleMutationalLoad mutationalLoad = new SampleMutationalLoad(snvBiallelic, snvNonBiallelic, indelBiallelic, indelNonBiallelic);

        return new SampleData(variants, mutationalLoad);
    }
}
