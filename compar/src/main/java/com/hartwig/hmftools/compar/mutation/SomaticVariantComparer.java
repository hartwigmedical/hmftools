package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.ComparConfig.NEW_SOURCE;
import static com.hartwig.hmftools.compar.ComparConfig.REF_SOURCE;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.INVALID_BOTH;
import static com.hartwig.hmftools.compar.MismatchType.INVALID_NEW;
import static com.hartwig.hmftools.compar.MismatchType.INVALID_REF;
import static com.hartwig.hmftools.compar.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.MismatchType.REF_ONLY;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_LPS;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.InvalidDataItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jooq.Record;
import org.jooq.Result;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

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
                variants.addAll(loadVariants(sourceSampleId, mConfig.DbConnections.get(sourceName)));
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

        if(!hasRefItems || !hasNewItems)
        {
            InvalidDataItem invalidDataItem = new InvalidDataItem(category());

            if(!hasRefItems && !hasNewItems)
                mismatches.add(new Mismatch(invalidDataItem, null, INVALID_BOTH, Collections.EMPTY_LIST));
            else if(!hasRefItems)
                mismatches.add(new Mismatch(invalidDataItem, null, INVALID_REF, Collections.EMPTY_LIST));
            else if(!hasNewItems)
                mismatches.add(new Mismatch(invalidDataItem, null, INVALID_NEW, Collections.EMPTY_LIST));

            return false;
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

                SomaticVariantData matchedVariant = null;
                MatchFilterStatus matchFilterStatus = null;

                int index2 = 0;
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
                    else if(newVariant.Position > refVariant.Position)
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
                        matchedVariant = unfilteredVariant;
                    }
                }

                if(matchedVariant != null)
                {
                    refVariants.remove(index1);

                    // skip checking for diffs if the items are not reportable
                    boolean eitherReportable = refVariant.reportable() || matchedVariant.reportable();

                    if(matchLevel != REPORTABLE || eitherReportable)
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

            refVariants.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                    .forEach(x -> mismatches.add(new Mismatch(x, null, REF_ONLY, emptyDiffs)));

            for(SomaticVariantData newVariant : newVariants)
            {
                if(matchLevel == REPORTABLE && !newVariant.reportable())
                    continue;

                SomaticVariantData unfilteredVariant = findUnfilteredVariant(newVariant, REF_SOURCE);

                if(unfilteredVariant != null)
                {
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
                    testVariant.Chromosome, testVariant.Position, ref, alt, VariantType.type(context),
                    "", false, Hotspot.fromVariant(context), VariantTier.fromContext(context),
                    false, "", "", "", "",
                    "", context.hasAttribute(LOCAL_PHASE_SET), (int)context.getPhredScaledQual(),
                    0, context.getFilters());
        }

        return null;
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

            if(variant.Gene.isEmpty())
                continue;

            if(mConfig.RestrictToDrivers && !mConfig.DriverGenes.contains(variant.Gene))
                continue;

            variants.add(variant);

            if(fileSources.RequiresLiftover && mConfig.LiftoverCache.hasMappings())
            {
                int newPosition = mConfig.LiftoverCache.convertPosition(variant.Chromosome, variant.Position);

                if(newPosition != UNMAPPED_POSITION)
                    variant.setComparisonCoordinates(RefGenomeFunctions.enforceChrPrefix(variant.Chromosome), newPosition);
            }
        }

        CMP_LOGGER.debug("sample({}) loaded {} somatic variants", sampleId, variants.size());

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
