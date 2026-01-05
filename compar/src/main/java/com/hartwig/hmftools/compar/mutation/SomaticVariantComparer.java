package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.compar.common.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.ComparConfig.NEW_SOURCE;
import static com.hartwig.hmftools.compar.ComparConfig.REF_SOURCE;
import static com.hartwig.hmftools.compar.common.CommonUtils.countsAsCalled;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_BOTH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_NEW;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_REF;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.REF_ONLY;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_LPS;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_PURITY_ADJUSTED_VAF;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TUMOR_SUPPORTING_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TUMOR_TOTAL_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_VARIANT_COPY_NUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.InvalidDataItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jooq.Record;
import org.jooq.Result;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariantComparer implements ItemComparer
{
    private final ComparConfig mConfig;
    private final Map<String,VcfFileReader> mUnfilteredVcfReaders;

    public SomaticVariantComparer(final ComparConfig config)
    {
        mConfig = config;
        mUnfilteredVcfReaders = Maps.newHashMap();
    }

    @Override
    public Category category() { return SOMATIC_VARIANT; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        // use a custom method optimised for large numbers of variants
        final MatchLevel matchLevel = mConfig.Categories.get(category());

        final List<SomaticVariantData> refVariants = loadVariants(sampleId, mConfig.SourceNames.get(0));
        final List<SomaticVariantData> newVariants = loadVariants(sampleId, mConfig.SourceNames.get(1));

        return identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel);
    }

    public boolean identifyMismatches(final String sampleId, final List<Mismatch> mismatches, final List<SomaticVariantData> refVariants,
            final List<SomaticVariantData> newVariants, final MatchLevel matchLevel)
    {
        boolean hasRefItems = refVariants != null;
        boolean hasNewItems = newVariants != null;
        final List<String> emptyDiffs = List.of();

        if(!hasRefItems || !hasNewItems)
        {
            InvalidDataItem invalidDataItem = new InvalidDataItem(category());

            if(!hasRefItems && !hasNewItems)
                mismatches.add(new Mismatch(invalidDataItem, null, INVALID_BOTH, emptyDiffs));
            else if(!hasRefItems)
                mismatches.add(new Mismatch(invalidDataItem, null, INVALID_REF, emptyDiffs));
            else if(!hasNewItems)
                mismatches.add(new Mismatch(invalidDataItem, null, INVALID_NEW, emptyDiffs));

            return false;
        }

        final String refSourceSampleId = mConfig.sourceSampleId(REF_SOURCE, sampleId);
        final String newSourceSampleId = mConfig.sourceSampleId(NEW_SOURCE, sampleId);

        final Map<String,List<SomaticVariantData>> refVariantsMap = buildVariantMap(refVariants);
        final Map<String,List<SomaticVariantData>> newVariantsMap = buildVariantMap(newVariants);
        final List<SomaticVariantData> emptyVariants = List.of();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = chromosome.toString();
            List<SomaticVariantData> chromosomeRefVariants = refVariantsMap.get(chrStr);
            List<SomaticVariantData> chromosomeNewVariants = newVariantsMap.get(chrStr);

            if(chromosomeNewVariants == null && chromosomeRefVariants == null)
                continue;

            if(chromosomeNewVariants == null)
                chromosomeNewVariants = emptyVariants;

            if(chromosomeRefVariants == null)
                chromosomeRefVariants = emptyVariants;

            int index1 = 0;
            int index2 = 0;
            while(index1 < chromosomeRefVariants.size())
            {
                final SomaticVariantData refVariant = chromosomeRefVariants.get(index1);

                SomaticVariantData matchedVariant = null;

                // shift index2 back to index at or before first potentially matching variant
                while(index2 > 0 && (index2 >= chromosomeNewVariants.size() || chromosomeNewVariants.get(index2).Position >= refVariant.comparisonPosition()))
                {
                    --index2;
                }

                while(index2 < chromosomeNewVariants.size())
                {
                    final SomaticVariantData newVariant = chromosomeNewVariants.get(index2);

                    if(refVariant.matches(newVariant))
                    {
                        matchedVariant = newVariant;
                        chromosomeNewVariants.remove(index2);
                        break;
                    }
                    else if(newVariant.Position > refVariant.comparisonPosition())
                    {
                        break;
                    }

                    ++index2;
                }

                if(matchedVariant == null)
                {
                    final SomaticVariantData unfilteredVariant = findUnfilteredVariant(refVariant, NEW_SOURCE, newSourceSampleId);

                    if(unfilteredVariant != null)
                    {
                        matchedVariant = unfilteredVariant;
                    }
                }

                if(matchedVariant != null)
                {
                    chromosomeRefVariants.remove(index1);

                    if(includeMismatchWithVariant(refVariant, matchLevel) || includeMismatchWithVariant(matchedVariant, matchLevel))
                    {
                        Mismatch mismatch = refVariant.findMismatch(matchedVariant, matchLevel, mConfig.Thresholds, mConfig.IncludeMatches);

                        if(mismatch != null)
                            mismatches.add(mismatch);
                    }
                }
                else
                {
                    ++index1;
                }
            }

            chromosomeRefVariants.stream().filter(x -> includeMismatchWithVariant(x, matchLevel))
                    .forEach(x -> mismatches.add(new Mismatch(x, null, REF_ONLY, emptyDiffs)));

            for(SomaticVariantData newVariant : chromosomeNewVariants)
            {
                if(!includeMismatchWithVariant(newVariant, matchLevel))
                    continue;

                SomaticVariantData unfilteredVariant = findUnfilteredVariant(newVariant, REF_SOURCE, refSourceSampleId);

                if(unfilteredVariant != null)
                {
                    mismatches.add(unfilteredVariant.findMismatch(newVariant, matchLevel, mConfig.Thresholds, mConfig.IncludeMatches));
                }
                else
                {
                    mismatches.add(new Mismatch(null, newVariant, NEW_ONLY, emptyDiffs));
                }
            }
        }

        return true;
    }

    protected SomaticVariantData findUnfilteredVariant(final SomaticVariantData testVariant, final String otherSource,
            final String sourceSampleId)
    {
        VcfFileReader unfilteredVcfReader = mUnfilteredVcfReaders.get(otherSource);

        if(unfilteredVcfReader == null)
            return null;

        List<VariantContext> candidates = unfilteredVcfReader.findVariants(
                testVariant.comparisonChromosome(), testVariant.comparisonPosition(), testVariant.comparisonPosition());

        for(VariantContext context : candidates)
        {
            String ref = context.getReference().getBaseString();
            String alt = !context.getAlternateAlleles().isEmpty() ? context.getAlternateAlleles().get(0).toString() : ref;

            if(!testVariant.Ref.equals(ref) || !testVariant.Alt.equals(alt))
                continue;

            final AllelicDepth tumorAllelicDepth = AllelicDepth.fromGenotype(context.getGenotype(sourceSampleId));
            return new SomaticVariantData(
                    context.getContig(), context.getStart(), ref, alt, VariantType.type(context),
                    "", false, Hotspot.fromVariant(context), VariantTier.fromContext(context),
                    false, "", "", "", "",
                    "", context.hasAttribute(LOCAL_PHASE_SET), (int)context.getPhredScaledQual(),
                    0, context.getFilters(), 0, 0,
                    tumorAllelicDepth.AlleleReadCount, tumorAllelicDepth.TotalReadCount, true, false,
                    testVariant.Chromosome, testVariant.Position);
        }

        return null;
    }

    private boolean includeMismatchWithVariant(SomaticVariantData variant, MatchLevel matchLevel)
    {
        return countsAsCalled(variant, matchLevel);
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
        thresholds.addFieldThreshold(FLD_VARIANT_COPY_NUMBER, 0.3, 0.15);
        thresholds.addFieldThreshold(FLD_PURITY_ADJUSTED_VAF, 0.2, 0);
        thresholds.addFieldThreshold(FLD_TUMOR_SUPPORTING_READ_COUNT, 1, 0.2);
        thresholds.addFieldThreshold(FLD_TUMOR_TOTAL_READ_COUNT, 1, 0.2);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        List<String> fieldNames = VariantCommon.comparedFieldNames();
        fieldNames.add(FLD_SUBCLONAL_LIKELIHOOD);
        fieldNames.add(FLD_LPS);
        return fieldNames;
    }

    private List<SomaticVariantData> loadVariants(final String sampleId, final String sourceName)
    {
        String sourceSampleId = mConfig.sourceSampleId(sourceName, sampleId);

        if(!mConfig.DbConnections.isEmpty())
        {
            return loadVariants(sourceSampleId, mConfig.DbConnections.get(sourceName), sourceName);
        }
        else
        {
            String sourceGermlineSampleId = mConfig.sourceGermlineSampleId(sourceName, sampleId);
            FileSources fileSources = mConfig.FileSources.get(sourceName);
            return loadVariants(sourceSampleId, FileSources.sampleInstance(fileSources, sourceSampleId, sourceGermlineSampleId));
        }
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        final List<ComparableItem> items = Lists.newArrayList();
        loadVariants(sampleId, dbAccess, sourceName).forEach(x -> items.add(x));
        return items;
    }

    private List<SomaticVariantData> loadVariants(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        final List<SomaticVariantData> variants = Lists.newArrayList();

        Result<Record> results = dbAccess.context()
                .select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .fetch();

        for(Record record : results)
        {
            final SomaticVariantData variant = SomaticVariantData.fromRecord(record, sourceName, mConfig);
            variants.add(variant);
        }

        return variants;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> items = Lists.newArrayList();
        loadVariants(sampleId, fileSources).forEach(x -> items.add(x));
        return items;
    }

    private List<SomaticVariantData> loadVariants(final String sampleId, final FileSources fileSources)
    {
        final List<SomaticVariantData> variants = Lists.newArrayList();

        // use the Purple suffix if not specified
        boolean usePurpleVcf = fileSources.SomaticVcf.isEmpty();
        String vcfFile = usePurpleVcf ? PurpleCommon.purpleSomaticVcfFile(fileSources.Purple, sampleId) : fileSources.SomaticVcf;

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        if(!vcfFileReader.fileValid())
        {
            CMP_LOGGER.error("failed to read somatic VCF file({})", vcfFile);
            return null;
        }

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(variantContext.isFiltered())
                continue;

            SomaticVariantData variant = SomaticVariantData.fromContext(variantContext, sampleId, false, usePurpleVcf, fileSources.Source, mConfig);

            if(mConfig.RestrictToDrivers && !mConfig.DriverGenes.contains(variant.Gene))
                continue;

            variants.add(variant);
        }

        CMP_LOGGER.debug("sample({}) loaded {} {} somatic variants", sampleId, fileSources.Source, variants.size());

        // prepare the unfiltered file source if configured
        if(!fileSources.SomaticUnfilteredVcf.isEmpty())
        {
            VcfFileReader unfilteredVcfReader = new VcfFileReader(fileSources.SomaticUnfilteredVcf);

            if(!unfilteredVcfReader.fileValid())
            {
                CMP_LOGGER.error("failed to read somatic unfiltered VCF file({})", fileSources.SomaticUnfilteredVcf);
                return null;
            }

            mUnfilteredVcfReaders.put(fileSources.Source, unfilteredVcfReader);
        }

        return variants;
    }
}
