package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CategoryType.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.countsAsCalled;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_BOTH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_NEW;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_OLD;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.OLD_ONLY;
import static com.hartwig.hmftools.compar.common.SourceType.NEW;
import static com.hartwig.hmftools.compar.common.SourceType.OLD;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_BIALLELIC;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_BIALLELIC_PROB;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_LPS;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.compar.mutation.VariantData.addComparerFields;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.HotspotType;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.InvalidDataItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceData;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariantComparer extends ItemComparer
{
    private final Map<SourceType,VcfFileReader> mUnfilteredVcfReaders;

    public SomaticVariantComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        addComparerFields(mFields, fieldCheckMap);

        mFields.add(new FieldInfo(
                FLD_BIALLELIC, getOrMakeFieldCheck(fieldCheckMap, FLD_BIALLELIC), null));

        mFields.add(new FieldInfo(
                FLD_BIALLELIC_PROB,
                getOrMakeFieldCheck(fieldCheckMap, FLD_BIALLELIC_PROB, 0.3, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                FLD_LPS, getOrMakeFieldCheck(fieldCheckMap, FLD_LPS), null));

        mFields.add(new FieldInfo(
                FLD_SUBCLONAL_LIKELIHOOD,
                getOrMakeFieldCheck(fieldCheckMap, FLD_SUBCLONAL_LIKELIHOOD, 0.6, null),
                "%.2f"));

        mUnfilteredVcfReaders = Maps.newHashMap();
    }

    @Override
    public CategoryType category() { return SOMATIC_VARIANT; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        // use a custom method optimised for large numbers of variants
        MatchLevel matchLevel = mConfig.MatchingLevel;

        List<SomaticVariantData> oldVariants = loadVariants(sampleId, OLD);
        List<SomaticVariantData> newVariants = loadVariants(sampleId, NEW);

        return identifyMismatches(sampleId, mismatches, oldVariants, newVariants, matchLevel);
    }

    public boolean identifyMismatches(
            final String sampleId, final List<Mismatch> mismatches, final List<SomaticVariantData> oldVariants,
            final List<SomaticVariantData> newVariants, final MatchLevel matchLevel)
    {
        boolean hasOldItems = oldVariants != null;
        boolean hasNewItems = newVariants != null;
        final List<String> emptyDiffs = List.of();

        if(!hasOldItems || !hasNewItems)
        {
            InvalidDataItem invalidDataItem = new InvalidDataItem(category());

            if(!hasOldItems && !hasNewItems)
                mismatches.add(new Mismatch(invalidDataItem, null, INVALID_BOTH, emptyDiffs));
            else if(!hasOldItems)
                mismatches.add(new Mismatch(invalidDataItem, null, INVALID_OLD, emptyDiffs));
            else if(!hasNewItems)
                mismatches.add(new Mismatch(invalidDataItem, null, INVALID_NEW, emptyDiffs));

            return false;
        }

        String oldSourceSampleId = mConfig.sourceSampleId(OLD, sampleId);
        String newSourceSampleId = mConfig.sourceSampleId(NEW, sampleId);

        final Map<String,List<SomaticVariantData>> oldVariantsMap = buildVariantMap(oldVariants);
        final Map<String,List<SomaticVariantData>> newVariantsMap = buildVariantMap(newVariants);
        final List<SomaticVariantData> emptyVariants = List.of();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = chromosome.toString();
            List<SomaticVariantData> chromosomeOldVariants = oldVariantsMap.get(chrStr);
            List<SomaticVariantData> chromosomeNewVariants = newVariantsMap.get(chrStr);

            if(chromosomeNewVariants == null && chromosomeOldVariants == null)
                continue;

            if(chromosomeNewVariants == null)
                chromosomeNewVariants = emptyVariants;

            if(chromosomeOldVariants == null)
                chromosomeOldVariants = emptyVariants;

            int index1 = 0;
            int index2 = 0;
            while(index1 < chromosomeOldVariants.size())
            {
                SomaticVariantData oldVariant = chromosomeOldVariants.get(index1);

                SomaticVariantData matchedVariant = null;

                // shift index2 back to index at or before first potentially matching variant
                while(index2 > 0 && (index2 >= chromosomeNewVariants.size() || chromosomeNewVariants.get(index2).Position >= oldVariant.comparisonPosition()))
                {
                    --index2;
                }

                while(index2 < chromosomeNewVariants.size())
                {
                    final SomaticVariantData newVariant = chromosomeNewVariants.get(index2);

                    if(oldVariant.matches(newVariant))
                    {
                        matchedVariant = newVariant;
                        chromosomeNewVariants.remove(index2);
                        break;
                    }
                    else if(newVariant.Position > oldVariant.comparisonPosition())
                    {
                        break;
                    }

                    ++index2;
                }

                if(matchedVariant == null)
                {
                    SomaticVariantData unfilteredVariant = findUnfilteredVariant(oldVariant, NEW, newSourceSampleId);

                    if(unfilteredVariant != null)
                    {
                        matchedVariant = unfilteredVariant;
                    }
                }

                if(matchedVariant != null)
                {
                    chromosomeOldVariants.remove(index1);

                    if(includeMismatchWithVariant(oldVariant, matchLevel) || includeMismatchWithVariant(matchedVariant, matchLevel))
                    {
                        Mismatch mismatch = oldVariant.findMismatch(this, matchedVariant, matchLevel, mConfig.IncludeMatches);

                        if(mismatch != null)
                            mismatches.add(mismatch);
                    }
                }
                else
                {
                    ++index1;
                }
            }

            chromosomeOldVariants.stream().filter(x -> includeMismatchWithVariant(x, matchLevel))
                    .forEach(x -> mismatches.add(new Mismatch(x, null, OLD_ONLY, emptyDiffs)));

            for(SomaticVariantData newVariant : chromosomeNewVariants)
            {
                if(!includeMismatchWithVariant(newVariant, matchLevel))
                    continue;

                SomaticVariantData unfilteredVariant = findUnfilteredVariant(newVariant, OLD, oldSourceSampleId);

                if(unfilteredVariant != null)
                {
                    mismatches.add(unfilteredVariant.findMismatch(this, newVariant, matchLevel, mConfig.IncludeMatches));
                }
                else
                {
                    mismatches.add(new Mismatch(null, newVariant, NEW_ONLY, emptyDiffs));
                }
            }
        }

        return true;
    }

    protected SomaticVariantData findUnfilteredVariant(
            final SomaticVariantData testVariant, final SourceType sourceType, final String sourceSampleId)
    {
        VcfFileReader unfilteredVcfReader = mUnfilteredVcfReaders.get(sourceType);

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

            // only extract fields set by Sage (ie not Pave nor Purple)
            return new SomaticVariantData(
                    context.getContig(), context.getStart(), ref, alt, VariantType.type(context),
                    "", false, HotspotType.fromVariant(context), VariantTier.fromContext(context),
                    "", "", "", "",
                    "", (int)context.getPhredScaledQual(), context.getFilters(), 0, 0,
                    tumorAllelicDepth.AlleleReadCount, tumorAllelicDepth.TotalReadCount, true,
                    testVariant.Chromosome, testVariant.Position,
                    false, 0, false, 0, mFields);
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
    public List<String> displayFieldNames()
    {
        List<String> fieldNames = VariantData.sharedFieldNames();
        fieldNames.add(FLD_BIALLELIC);
        fieldNames.add(FLD_BIALLELIC_PROB);
        fieldNames.add(FLD_LPS);
        fieldNames.add(FLD_SUBCLONAL_LIKELIHOOD);
        return fieldNames;
    }

    private List<SomaticVariantData> loadVariants(final String sampleId, final SourceType sourceType)
    {
        String sourceSampleId = mConfig.sourceSampleId(sourceType, sampleId);
        SourceData sourceData = mConfig.getSourceData(sourceType);

        String sourceReferenceId = mConfig.sourceReferenceId(sourceType, sampleId);
        return loadVariants(sourceSampleId, PipelineSourcePaths.sampleInstance(sourceData.PipelinePaths, sourceSampleId, sourceReferenceId));
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> items = Lists.newArrayList();
        loadVariants(sampleId, fileSources).forEach(x -> items.add(x));
        return items;
    }

    private List<SomaticVariantData> loadVariants(final String sampleId, final PipelineSourcePaths fileSources)
    {
        final List<SomaticVariantData> variants = Lists.newArrayList();

        String vcfFile = PurpleCommon.purpleSomaticVcfFile(fileSources.Purple, sampleId);

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        if(!vcfFileReader.fileValid())
        {
            CMP_LOGGER.warn("failed to read somatic VCF file({})", vcfFile);
            return null;
        }

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(variantContext.isFiltered())
                continue;

            SomaticVariantData variant = SomaticVariantData.fromContext(
                    variantContext, sampleId, false, fileSources.Source, mConfig, mFields);

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
