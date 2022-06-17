package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.MismatchType.REF_ONLY;
import static com.hartwig.hmftools.compar.somatic.SomaticVariantData.FLD_LPS;
import static com.hartwig.hmftools.compar.somatic.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.compar.somatic.VariantCommon.FLD_QUAL;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

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
        // use a custom method optimised for large numbers of variants
        final MatchLevel matchLevel = mConfig.Categories.get(category());

        final List<SomaticVariantData> allRefVariants = Lists.newArrayList();
        final List<SomaticVariantData> allNewVariants = Lists.newArrayList();

        for(int i = 0; i <= 1; ++i)
        {
            final List<SomaticVariantData> variants = (i == 0) ? allRefVariants : allNewVariants;

            final String sourceName = mConfig.SourceNames.get(i);
            String sourceSampleId = mConfig.sourceSampleId(sourceName, sampleId);

            if(!mConfig.DbConnections.isEmpty())
            {
                variants.addAll(loadVariants(sourceSampleId, mConfig.DbConnections.get(sourceName)));
            }
            else
            {
                FileSources fileSources = mConfig.FileSources.get(sourceName);
                variants.addAll(loadVariants(sourceSampleId, FileSources.sampleInstance(fileSources, sampleId)));
            }
        }

        final Map<String,List<SomaticVariantData>> refVariantsMap = buildVariantMap(allRefVariants);
        final Map<String,List<SomaticVariantData>> newVariantsMap = buildVariantMap(allNewVariants);

        final List<String> emptyDiffs = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = chromosome.toString();
            List<SomaticVariantData> refVariants = refVariantsMap.get(chrStr);
            List<SomaticVariantData> newVariants = newVariantsMap.get(chrStr);

            if(newVariants == null && refVariants == null)
                continue;

            if(newVariants == null && refVariants != null)
            {
                refVariants.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                        .forEach(x -> mismatches.add(new Mismatch(x, null, REF_ONLY, emptyDiffs)));
                continue;
            }
            else if(refVariants == null && newVariants != null)
            {
                newVariants.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                        .forEach(x -> mismatches.add(new Mismatch(null, x, NEW_ONLY, emptyDiffs)));
                continue;
            }

            int index1 = 0;
            while(index1 < refVariants.size())
            {
                final SomaticVariantData refVariant = refVariants.get(index1);

                boolean matched = false;

                int index2 = 0;
                while(index2 < newVariants.size())
                {
                    final SomaticVariantData newVariant = newVariants.get(index2);

                    if(refVariant.matches(newVariant))
                    {
                        refVariants.remove(index1);
                        newVariants.remove(index2);
                        matched = true;

                        // skip checking for diffs if the items are not reportable
                        boolean eitherReportable = newVariant.reportable() || newVariant.reportable();

                        if(matchLevel != REPORTABLE || eitherReportable)
                        {
                            Mismatch mismatch = refVariant.findMismatch(newVariant, matchLevel, mConfig.Thresholds);

                            if(mismatch != null)
                                mismatches.add(mismatch);
                        }

                        break;
                    }
                    else
                    {
                        ++index2;
                    }
                }

                if(!matched)
                    ++index1;
            }

            refVariants.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                    .forEach(x -> mismatches.add(new Mismatch(x, null, REF_ONLY, emptyDiffs)));

            newVariants.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                    .forEach(x -> mismatches.add(new Mismatch(null, x, NEW_ONLY, emptyDiffs)));
        }
    }

    private Map<String,List<SomaticVariantData>> buildVariantMap(final List<SomaticVariantData> variants)
    {
        final Map<String,List<SomaticVariantData>> chrVariantsMap = Maps.newHashMap();

        for(SomaticVariantData variant : variants)
        {
            String chromosome = RefGenomeFunctions.stripChrPrefix(variant.Chromosome);
            List<SomaticVariantData> chrVariants = chrVariantsMap.get(chromosome);
            if(chrVariants == null)
            {
                chrVariants = Lists.newArrayList();
                chrVariantsMap.put(chromosome, chrVariants);
            }

            chrVariants.add(variant);
        }

        return chrVariantsMap;
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_QUAL, 20, 0.2);
        thresholds.addFieldThreshold(FLD_SUBCLONAL_LIKELIHOOD, 0.6, 0);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        List<String> fieldNames = VariantCommon.comparedFieldNames();
        fieldNames.add(FLD_SUBCLONAL_LIKELIHOOD);
        fieldNames.add(FLD_LPS);
        return fieldNames;
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<ComparableItem> items = Lists.newArrayList();
        loadVariants(sampleId, dbAccess).forEach(x -> items.add(x));
        return items;
    }

    private List<SomaticVariantData> loadVariants(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<SomaticVariantData> variants = Lists.newArrayList();

        Result<Record> results = dbAccess.context()
                .select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.GENE.isNotNull())
                .fetch();

        for(Record record : results)
        {
            variants.add(SomaticVariantData.fromRecord(record));
        }

        return variants;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> items = Lists.newArrayList();
        loadVariants(sampleId, fileSources).forEach(x -> items.add(x));
        return items;
    }

    private List<SomaticVariantData> loadVariants(final String sampleId, final FileSources fileSources)
    {
        final List<SomaticVariantData> variants = Lists.newArrayList();

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        // use the Purple suffix if not specified
        String vcfFile = !fileSources.SomaticVcf.isEmpty() ?
                fileSources.SomaticVcf : fileSources.Purple + sampleId + PURPLE_SOMATIC_VCF_SUFFIX;

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);

            for(VariantContext variantContext : reader.iterator())
            {
                if(!filter.test(variantContext))
                    continue;

                SomaticVariantData variant = SomaticVariantData.fromContext(variantContext);

                if(!variant.Gene.isEmpty())
                    variants.add(variant);
            }

            CMP_LOGGER.debug("sample({}) loaded {} somatic variants", sampleId, variants.size());
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to read somatic VCF file({}): {}", vcfFile, e.toString());
        }

        return variants;
    }
}
