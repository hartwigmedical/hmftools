package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.sage.SageCommon.PAVE_SOMATIC_VCF_ID;
import static com.hartwig.hmftools.common.sage.SageCommon.SAGE_FILE_ID;
import static com.hartwig.hmftools.common.sage.SageCommon.SAGE_SOMATIC_VCF_ID;
import static com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser.fromVariantContext;
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

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
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
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariantComparer extends ItemComparer
{
    private final Map<SourceType,VcfFileReader> mSageVcfReaders;
    private final Map<SourceType,VcfFileReader> mPaveVcfReaders;

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

        mSageVcfReaders = Maps.newHashMap();
        mPaveVcfReaders = Maps.newHashMap();
    }

    @Override
    public CategoryType category() { return SOMATIC_VARIANT; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        mSageVcfReaders.clear();
        mPaveVcfReaders.clear();

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
        List<String> emptyDiffs = List.of();

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

        boolean oldTruthsetSourced = mConfig.isTruthsetSourced(SourceType.OLD);
        boolean newTruthsetSourced = mConfig.isTruthsetSourced(SourceType.NEW);

        Map<String,List<SomaticVariantData>> oldVariantsMap = buildVariantMap(oldVariants);
        Map<String,List<SomaticVariantData>> newVariantsMap = buildVariantMap(newVariants);
        List<SomaticVariantData> emptyVariants = List.of();

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
                while(index2 > 0
                && (index2 >= chromosomeNewVariants.size() || chromosomeNewVariants.get(index2).Position >= oldVariant.Position))
                {
                    --index2;
                }

                while(index2 < chromosomeNewVariants.size())
                {
                    SomaticVariantData newVariant = chromosomeNewVariants.get(index2);

                    if(oldVariant.matches(newVariant))
                    {
                        matchedVariant = newVariant;
                        chromosomeNewVariants.remove(index2);
                        break;
                    }
                    else if(newVariant.Position > oldVariant.Position)
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

                    if((includeMismatchWithVariant(oldVariant, matchLevel) || oldTruthsetSourced)
                    || (includeMismatchWithVariant(matchedVariant, matchLevel) || newTruthsetSourced))
                    {
                        Mismatch mismatch = oldVariant.findMismatch(
                                this, matchedVariant, matchLevel, mConfig.IncludeMatches,
                                oldTruthsetSourced || newTruthsetSourced);

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
                if(!newTruthsetSourced && !includeMismatchWithVariant(newVariant, matchLevel))
                    continue;

                SomaticVariantData unfilteredVariant = findUnfilteredVariant(newVariant, OLD, oldSourceSampleId);

                if(unfilteredVariant != null)
                {
                    mismatches.add(unfilteredVariant.findMismatch(
                            this, newVariant, matchLevel, mConfig.IncludeMatches, newTruthsetSourced));
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
        // first look in the Pave VCFs, then back into the Sage ones
        for(int i = 0; i <= 1; ++i)
        {
            VcfFileReader vcfReader = (i == 0) ? mPaveVcfReaders.get(sourceType) : mSageVcfReaders.get(sourceType);

            if(vcfReader == null)
                continue;

            List<VariantContext> candidates = vcfReader.findVariants(
                    testVariant.Chromosome, testVariant.Position, testVariant.Position);

            for(VariantContext variantContext : candidates)
            {
                String ref = variantContext.getReference().getBaseString();
                String alt = !variantContext.getAlternateAlleles().isEmpty() ? variantContext.getAlternateAlleles().get(0).toString() : ref;

                if(!testVariant.Ref.equals(ref) || !testVariant.Alt.equals(alt))
                    continue;

                return SomaticVariantData.fromContext(variantContext, sourceSampleId, true, sourceType, mConfig, mFields);
            }
        }

        return null;
    }

    private boolean includeMismatchWithVariant(final SomaticVariantData variant, final MatchLevel matchLevel)
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

        List<SomaticVariantData> variants = Lists.newArrayList();

        if(sourceData.PipelinePaths != null)
        {
            List<SomaticVariantData> fileVariants = loadVariants(
                    sourceSampleId, PipelineSourcePaths.sampleInstance(sourceData.PipelinePaths, sourceSampleId, sourceReferenceId));
            if(fileVariants == null)
            {
                return null;
            }
            variants.addAll(fileVariants);
        }
        else
        {
            List<TruthsetValue> truthsetValues = sourceData.Truthset.sampleTruthsetEntries(sampleId, category());

            if(truthsetValues == null || truthsetValues.isEmpty())
                return Collections.emptyList();

            // validate truthset fields against the comparer
            List<FieldInfo> fields = fieldsList();

            for(TruthsetValue truthsetValue : truthsetValues)
            {
                if(fields.stream().noneMatch(x -> x.Name.equals(truthsetValue.FieldName)))
                {
                    CMP_LOGGER.error("category({}) invalid truthset entry({})", category(), truthsetValue);
                    return null;
                }
            }

            // group by keys
            Map<String, List<TruthsetValue>> truthsetValuesByKey = truthsetValues.stream()
                    .collect(Collectors.groupingBy(x -> x.Key));

            for(Map.Entry<String,List<TruthsetValue>> entry : truthsetValuesByKey.entrySet())
            {
                SomaticVariantData variant = SomaticVariantData.fromTruthset(entry.getValue(), mFields);

                if(variant == null)
                    return null;

                variants.add(variant);
            }
        }

        return variants;
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
        List<SomaticVariantData> variants = Lists.newArrayList();

        String vcfFile = PurpleCommon.purpleSomaticVcfFile(fileSources.Purple, sampleId);

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        if(!vcfFileReader.fileValid())
        {
            CMP_LOGGER.warn("failed to read somatic VCF file({})", vcfFile);
            return null;
        }

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
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
            boolean isSage = fileSources.SomaticUnfilteredVcf.contains(SAGE_FILE_ID);

            loadVcfReader(
                    fileSources.SomaticUnfilteredVcf,
                    isSage ? mSageVcfReaders : mPaveVcfReaders, fileSources.Source);
        }
        else
        {
            String sageVcf = fileSources.SageSomatic + sampleId + SAGE_SOMATIC_VCF_ID;
            String paveVcf = fileSources.PaveSomatic + sampleId + PAVE_SOMATIC_VCF_ID;

            if(!fileSources.SageSomatic.isEmpty() && Files.exists(Paths.get(sageVcf)))
            {
                loadVcfReader(sageVcf, mSageVcfReaders, fileSources.Source);
            }

            if(!fileSources.PaveSomatic.isEmpty() && Files.exists(Paths.get(paveVcf)))
            {
                loadVcfReader(paveVcf, mPaveVcfReaders, fileSources.Source);
            }
        }

        return variants;
    }

    private static void loadVcfReader(final String vcfFilename, final Map<SourceType,VcfFileReader> vcfReaders, final SourceType source)
    {
        VcfFileReader vcfReader = new VcfFileReader(vcfFilename);

        if(!vcfReader.fileValid())
        {
            CMP_LOGGER.error("failed to read somatic unfiltered VCF file({})", vcfFilename);
            return;
        }

        vcfReaders.put(source, vcfReader);
    }
}
