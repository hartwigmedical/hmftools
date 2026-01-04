package com.hartwig.hmftools.sage.phase;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LPS_APPEND_INFO;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class AppendVariantPhaser implements VariantPhaser
{
    private final Map<Integer,LocalPhaseSet> mLocalPhaseSets;
    private final Map<SimpleVariant,List<Integer>> mVariantLpsIds;
    private final Map<String,Map<Integer,LpsReadCounts>> mSampleLpsReadCounts;

    private Map<Integer, LpsReadCounts> mCurrentSampleLpsCounts;

    public AppendVariantPhaser()
    {
        mLocalPhaseSets = Maps.newHashMap();
        mSampleLpsReadCounts = Maps.newHashMap();
        mVariantLpsIds = Maps.newHashMap();
        mCurrentSampleLpsCounts = null;
    }

    @Override
    public void initialise(final ChrBaseRegion region, final String sample)
    {
        mCurrentSampleLpsCounts = Maps.newHashMap();
        mSampleLpsReadCounts.put(sample, mCurrentSampleLpsCounts);
    }

    @Override
    public void registeredPhasedVariants(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        if(posCounters.size() >= 2)
        {
            List<SimpleVariant> supportVariants = posCounters.stream().map(x -> x.variant()).collect(Collectors.toList());

            if(!supportVariants.isEmpty())
                registeredPhasedVariants(supportVariants, true);
        }

        if(negCounters.size() >= 2)
        {
            List<SimpleVariant> depthVariants = negCounters.stream().map(x -> x.variant()).collect(Collectors.toList());

            if(!depthVariants.isEmpty())
                registeredPhasedVariants(depthVariants, false);
        }
    }

    private void registeredPhasedVariants(final List<SimpleVariant> variants, boolean hasAltSupport)
    {
        Set<Integer> processedLpsIds = Sets.newHashSet();

        for(SimpleVariant variant : variants)
        {
            List<Integer> lpsIds = mVariantLpsIds.get(variant);

            if(lpsIds == null)
                continue;

            for(Integer lpsId : lpsIds)
            {
                if(processedLpsIds.contains(lpsId))
                    continue;

                processedLpsIds.add(lpsId);

                LocalPhaseSet localPhaseSet = mLocalPhaseSets.get(lpsId);

                if(localPhaseSet == null)
                {
                    SG_LOGGER.error("lpsId({}) missing", lpsId);
                    continue;
                }

                List<SimpleVariant> lpsVariants = variants;

                // skip if the variants from the read aren't an exact match for this LPS and are for another
                if(!localPhaseSet.variantsMatch(variants) && variantsMatchLps(variants))
                    continue;

                if(variants.stream().anyMatch(x -> !localPhaseSet.Variants.contains(x)))
                {
                    // remove any variants not in the LPS
                    lpsVariants = variants.stream().filter(x -> localPhaseSet.Variants.contains(x)).collect(Collectors.toList());

                    if(lpsVariants.size() < 2)
                        continue;
                }

                LpsReadCounts readCounts = mCurrentSampleLpsCounts.get(lpsId);

                if(readCounts == null)
                {
                    readCounts = new LpsReadCounts(lpsId, localPhaseSet.Variants);
                    mCurrentSampleLpsCounts.put(lpsId, readCounts);
                }


                // only look for subsets if this variant set doesn't match any others precisely
                if(lpsVariants.size() < readCounts.Variants.size() && variantsMatchLps(lpsVariants))
                    continue;

                readCounts.registerCount(lpsVariants, hasAltSupport);
            }
        }
    }

    private boolean variantsMatchLps(final List<SimpleVariant> variants)
    {
        for(LocalPhaseSet localPhaseSet : mLocalPhaseSets.values())
        {
            if(localPhaseSet.variantsMatch(variants))
                return true;
        }

        return false;
    }

    private class LocalPhaseSet
    {
        public final int Id;
        public final List<SimpleVariant> Variants;

        public LocalPhaseSet(final int id, final List<SimpleVariant> variants)
        {
            Id = id;
            Variants = variants;
        }

        public int variantCount() { return Variants.size(); }

        public boolean variantsMatch(final List<SimpleVariant> variants)
        {
            if(Variants.size() != variants.size())
                return false;

            return Variants.stream().allMatch(x -> variants.contains(x));
        }

        public String toString()
        {
            return format("%d: variants(%d)", Id, Variants.size());
        }
    }

    public void registerLocalPhaseSets(final List<Candidate> candidates, final List<VariantContext> variantContexts)
    {
        for(int i = 0; i < candidates.size(); ++i)
        {
            VariantContext variantContext = variantContexts.get(i);

            if(variantContext.isFiltered() && !variantContext.getFilters().contains(PASS_FILTER))
                continue;

            if(!variantContext.hasAttribute(LOCAL_PHASE_SET))
                continue;

            Candidate candidate = candidates.get(i);

            List<Integer> lpsIds = variantContext.getAttributeAsIntList(LOCAL_PHASE_SET, 0);

            for(int lpsId : lpsIds)
            {
                LocalPhaseSet localPhaseSet = mLocalPhaseSets.get(lpsId);

                if(localPhaseSet == null)
                {
                    localPhaseSet = new LocalPhaseSet(lpsId, Lists.newArrayList(candidate.variant()));
                    mLocalPhaseSets.put(lpsId, localPhaseSet);
                }
                else
                {
                    localPhaseSet.Variants.add(candidate.variant());
                }
            }
        }

        // purge any single-variant groups
        Set<Integer> singleVariantIds = mLocalPhaseSets.entrySet().stream()
                .filter(x -> x.getValue().variantCount() == 1).map(x -> x.getKey()).collect(Collectors.toSet());

        singleVariantIds.forEach(x -> mLocalPhaseSets.remove(x));

        // make links from variants to their LPS IDs
        for(LocalPhaseSet localPhaseSet : mLocalPhaseSets.values())
        {
            for(SimpleVariant variant : localPhaseSet.Variants)
            {
                List<Integer> variantLpsIds = mVariantLpsIds.get(variant);

                if(variantLpsIds == null)
                {
                    variantLpsIds = Lists.newArrayList();
                    mVariantLpsIds.put(variant, variantLpsIds);
                }

                variantLpsIds.add(localPhaseSet.Id);
            }
        }
    }

    public void populateLocalPhaseSetInfo(final List<Candidate> candidates, final List<VariantContext> variantContexts)
    {
        boolean logAltOnly = false;

        for(int i = 0; i < candidates.size(); ++i)
        {
            Candidate candidate = candidates.get(i);

            List<Integer> lpsIds = mVariantLpsIds.get(candidate.variant());

            if(lpsIds == null)
                continue;

            VariantContext variantContext = variantContexts.get(i);

            for(Genotype genotype : variantContext.getGenotypes())
            {
                Map<Integer,LpsReadCounts> sampleReadCounts = mSampleLpsReadCounts.get(genotype.getSampleName());

                if(sampleReadCounts == null)
                    continue;

                StringJoiner lpsInfo = new StringJoiner(ITEM_DELIM);

                for(int lpsId : lpsIds)
                {
                    LpsReadCounts readCounts = sampleReadCounts.get(lpsId);

                    if(readCounts == null)
                        continue;

                    if(!logAltOnly || readCounts.AltSupport > 0) // during testing only log if has alt support
                    {
                        lpsInfo.add(format("%d=%d/%d", readCounts.Id, readCounts.AltSupport, readCounts.Depth));
                    }

                    if(readCounts.SubsetVariantCounts != null)
                    {
                        for(LpsReadCounts subsetReadCounts : readCounts.SubsetVariantCounts)
                        {
                            if(subsetReadCounts.Variants.contains(candidate.variant())
                            && (!logAltOnly || subsetReadCounts.AltSupport > 0))
                            {
                                lpsInfo.add(format("%d_%d=%d/%d",
                                        readCounts.Id, subsetReadCounts.Id, subsetReadCounts.AltSupport, subsetReadCounts.Depth));
                            }
                        }
                    }
                }

                genotype.getExtendedAttributes().put(LPS_APPEND_INFO, lpsInfo.toString());
            }
        }
    }

    @VisibleForTesting
    public Map<String,Map<Integer,LpsReadCounts>> sampleLpsReadCounts() { return mSampleLpsReadCounts; }
}