package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.compar.common.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.ComparConfig.NEW_SOURCE;
import static com.hartwig.hmftools.compar.ComparConfig.REF_SOURCE;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_BOTH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_NEW;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_REF;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.REF_ONLY;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_LPS;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.region.BasePosition;
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

        final List<SomaticVariantData> allRefVariants = Lists.newArrayList();
        final List<SomaticVariantData> allNewVariants = Lists.newArrayList();

        boolean usesNonPurpleVcfs = false;
        boolean hasRefItems = false;
        boolean hasNewItems = false;

        for(int i = 0; i <= 1; ++i)
        {
            final List<SomaticVariantData> variants = (i == 0) ? allRefVariants : allNewVariants;

            final String sourceName = mConfig.SourceNames.get(i);

            String sourceSampleId = mConfig.sourceSampleId(sourceName, sampleId);

            if(!mConfig.DbConnections.isEmpty())
            {
                variants.addAll(loadVariants(sourceSampleId, mConfig.DbConnections.get(sourceName), sourceName));
            }
            else
            {
                FileSources fileSources = mConfig.FileSources.get(sourceName);
                List<SomaticVariantData> fileVariants = loadVariants(sourceSampleId, FileSources.sampleInstance(fileSources, sourceSampleId));

                if(fileVariants == null)
                    continue;

                variants.addAll(fileVariants);
                usesNonPurpleVcfs |= !fileSources.SomaticVcf.isEmpty();
            }

            if(sourceName.equals(REF_SOURCE))
                hasRefItems = true;
            else
                hasNewItems = true;
        }

        final List<String> emptyDiffs = Lists.newArrayList();
        
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

        final Map<String,List<SomaticVariantData>> refVariantsMap = buildVariantMap(allRefVariants);
        final Map<String,List<SomaticVariantData>> newVariantsMap = buildVariantMap(allNewVariants);
        final List<SomaticVariantData> emptyVariants = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = chromosome.toString();
            List<SomaticVariantData> refVariants = refVariantsMap.get(chrStr);
            List<SomaticVariantData> newVariants = newVariantsMap.get(chrStr);

            if(newVariants == null && refVariants == null)
                continue;

            if(newVariants == null)
                newVariants = emptyVariants;

            if(refVariants == null)
                refVariants = emptyVariants;

            int index1 = 0;
            int index2 = 0;
            while(index1 < refVariants.size())
            {
                final SomaticVariantData refVariant = refVariants.get(index1);

                SomaticVariantData matchedVariant = null;
                MatchFilterStatus matchFilterStatus = null;

                // shift index2 back to index at or before first potentially matching variant
                while(index2 > 0 && (index2 >= newVariants.size() || newVariants.get(index2).Position >= refVariant.comparisonPosition()))
                {
                    --index2;
                }

                while(index2 < newVariants.size())
                {
                    final SomaticVariantData newVariant = newVariants.get(index2);

                    if(refVariant.matches(newVariant))
                    {
                        matchedVariant = newVariant;
                        matchFilterStatus = MatchFilterStatus.BOTH_UNFILTERED;
                        newVariants.remove(index2);
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
                    final SomaticVariantData unfilteredVariant = findUnfilteredVariant(refVariant, NEW_SOURCE);

                    if(unfilteredVariant != null)
                    {
                        matchFilterStatus = MatchFilterStatus.NEW_FILTERED;
                        unfilteredVariant.setComparisonCoordinates(refVariant.Chromosome, refVariant.Position);
                        matchedVariant = unfilteredVariant;
                    }
                }

                if(matchedVariant != null)
                {
                    refVariants.remove(index1);

                    if(includeMismatchWithVariant(refVariant, matchLevel) || includeMismatchWithVariant(matchedVariant, matchLevel))
                    {
                        Mismatch mismatch = refVariant.findDiffs(matchedVariant, mConfig.Thresholds, matchFilterStatus, usesNonPurpleVcfs);

                        if(mismatch != null)
                            mismatches.add(mismatch);
                    }
                }
                else
                {
                    ++index1;
                }
            }

            refVariants.stream().filter(x -> includeMismatchWithVariant(x, matchLevel))
                    .forEach(x -> mismatches.add(new Mismatch(x, null, REF_ONLY, emptyDiffs)));

            for(SomaticVariantData newVariant : newVariants)
            {
                if(!includeMismatchWithVariant(newVariant, matchLevel))
                    continue;

                SomaticVariantData unfilteredVariant = findUnfilteredVariant(newVariant, REF_SOURCE);

                if(unfilteredVariant != null)
                {
                    unfilteredVariant.setComparisonCoordinates(newVariant.Chromosome, newVariant.Position);
                    mismatches.add(unfilteredVariant.findDiffs(newVariant, mConfig.Thresholds, MatchFilterStatus.REF_FILTERED, usesNonPurpleVcfs));
                }
                else
                {
                    mismatches.add(new Mismatch(null, newVariant, NEW_ONLY, emptyDiffs));
                }
            }
        }

        return true;
    }

    protected SomaticVariantData findUnfilteredVariant(final SomaticVariantData testVariant, final String otherSource)
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

            return new SomaticVariantData(
                    context.getContig(), context.getStart(), ref, alt, VariantType.type(context),
                    "", false, Hotspot.fromVariant(context), VariantTier.fromContext(context),
                    false, "", "", "", "",
                    "", context.hasAttribute(LOCAL_PHASE_SET), (int)context.getPhredScaledQual(),
                    0, context.getFilters());
        }

        return null;
    }

    private boolean includeMismatchWithVariant(SomaticVariantData variant, MatchLevel matchLevel)
    {
        boolean reportabilityIsFine = (matchLevel != REPORTABLE || variant.reportable());
        boolean isInGene = !variant.Gene.isEmpty();
        return reportabilityIsFine && isInGene;
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
            final SomaticVariantData variant = SomaticVariantData.fromRecord(record);
            BasePosition comparisonPosition = determineComparisonGenomePosition(
                    variant.Chromosome, variant.Position, sourceName, mConfig.RequiresLiftover, mConfig.LiftoverCache);
            variant.setComparisonCoordinates(comparisonPosition.Chromosome, comparisonPosition.Position);
            variants.add(variant);
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

        // use the Purple suffix if not specified
        String vcfFile = !fileSources.SomaticVcf.isEmpty() ?
                fileSources.SomaticVcf : PurpleCommon.purpleSomaticVcfFile(fileSources.Purple, sampleId);

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

            SomaticVariantData variant = SomaticVariantData.fromContext(variantContext);

            if(mConfig.RestrictToDrivers && !mConfig.DriverGenes.contains(variant.Gene))
                continue;

            BasePosition comparisonPosition = determineComparisonGenomePosition(
                    variant.Chromosome, variant.Position, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);
            variant.setComparisonCoordinates(comparisonPosition.Chromosome, comparisonPosition.Position);

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
