package com.hartwig.hmftools.dnds.builder;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;
import static com.hartwig.hmftools.dnds.DndsCommon.MAX_REPEAT_COUNT;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.byteToBoolean;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import static org.jooq.impl.DSL.count;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.region.BedFileReader;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
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
    private final Map<String, List<ChrBaseRegion>> mTargetRegionsByChromosome;

    public SampleDataLoader(final String purpleDir, final DatabaseAccess dbAccess, final String targetRegionsBedFile)
    {
        mDbAccess = dbAccess;
        mPurpleDir = purpleDir;
        mTargetRegionsByChromosome = loadTargetRegionsByChromosome(targetRegionsBedFile);
    }

    public SampleData loadSampleData(final String sampleId)
    {
        if(mDbAccess != null)
        {
            return new SampleData(processVariantsFromDatabase(sampleId, MAX_REPEAT_COUNT), readMutationalLoadFromDatabase(sampleId));
        }
        else if(mPurpleDir != null)
        {
            return processVariantsFromVcf(sampleId, MAX_REPEAT_COUNT);
        }
        else
        {
            DN_LOGGER.error("sample({}) neither database nor purple directory configured", sampleId);
            return new SampleData(Collections.emptyList(), null);
        }
    }

    private static Map<String, List<ChrBaseRegion>> loadTargetRegionsByChromosome(final String bedFile)
    {
        DN_LOGGER.info("loading bed file: {}", bedFile);

        List<ChrBaseRegion> targetRegions;

        try
        {
            targetRegions = BedFileReader.loadBedFile(bedFile);
        }
        catch(Exception e)
        {
            DN_LOGGER.error("failed to load bed file({}): {}", bedFile, e.toString());
            System.exit(1);
            return null;
        }

        Map<String, List<ChrBaseRegion>> regionsByChromosome = Maps.newHashMap();

        for(ChrBaseRegion region : targetRegions)
        {
            regionsByChromosome.computeIfAbsent(region.Chromosome, k -> Lists.newArrayList()).add(region);
        }

        return regionsByChromosome;
    }

    @VisibleForTesting
    static boolean inTargetRegions(final Map<String, List<ChrBaseRegion>> targetRegionsByChromosome, final String chromosome, final int position)
    {
        // Use binary search to speed up VCF slicing - target regions are sorted by chromosome then start position

        // normalise since target regions are keyed without a 'chr' prefix (see BedLine.toRegion), whereas VCF contigs
        // may or may not carry one depending on ref genome version (see RefGenomeVersion.versionedChromosome)
        if(!HumanChromosome.contains(chromosome))
            return false;

        List<ChrBaseRegion> regions = targetRegionsByChromosome.get(HumanChromosome.fromString(chromosome).toString());

        if(regions == null)
            return false;

        int lowerIndex = 0;
        int upperIndex = regions.size() - 1;

        while(lowerIndex <= upperIndex)
        {
            int midIndex = (lowerIndex + upperIndex) / 2;
            ChrBaseRegion region = regions.get(midIndex);

            if(position < region.start())
                upperIndex = midIndex - 1;
            else if(position > region.end())
                lowerIndex = midIndex + 1;
            else
                return true;
        }

        return false;
    }

    private SampleMutationalLoad readMutationalLoadFromDatabase(final String sample)
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

    private List<SomaticVariant> processVariantsFromDatabase(final String sample, int maxRepeatCount)
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
            String chromosome = record.getValue(SOMATICVARIANT.CHROMOSOME);
            int position = record.getValue(SOMATICVARIANT.POSITION);

            if(!inTargetRegions(mTargetRegionsByChromosome, chromosome, position))
                continue;

            SomaticVariant variant = new SomaticVariant(
                    chromosome,
                    position,
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

    private SampleData processVariantsFromVcf(final String sample, final int maxRepeatCount)
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

                int repeatCount = decorator.repeatCount();
                if(decorator.repeatCount() > maxRepeatCount)
                    continue;

                if(!inTargetRegions(mTargetRegionsByChromosome, decorator.chromosome(), decorator.position()))
                    continue;

                VariantImpact impact = decorator.variantImpact();
                SomaticVariant variant = new SomaticVariant(
                        decorator.chromosome(),
                        decorator.position(),
                        decorator.ref(),
                        decorator.alt(),
                        decorator.gene(),
                        isBiallelic,
                        decorator.isHotspot(),
                        impact.WorstCodingEffect,
                        impact.CanonicalCodingEffect,
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
