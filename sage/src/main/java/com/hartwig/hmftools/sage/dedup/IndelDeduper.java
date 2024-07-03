package com.hartwig.hmftools.sage.dedup;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.INDEL_DEDUP_MIN_MATCHED_LPS_PERCENT;
import static com.hartwig.hmftools.sage.SageConstants.MAX_READ_EDGE_DISTANCE_PERC;
import static com.hartwig.hmftools.sage.filter.SoftFilter.DEDUP_INDEL;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.filter.SoftFilter;

public class IndelDeduper
{
    private final RefGenomeInterface mRefGenome;
    private int mGroupIterations;
    private final int mReadEdgeDistanceThreshold;

    public IndelDeduper(final RefGenomeInterface refGenome, int readLength)
    {
        mRefGenome = refGenome;
        mReadEdgeDistanceThreshold = (int)(readLength * MAX_READ_EDGE_DISTANCE_PERC);
        mGroupIterations = 0;
    }

    public void dedupVariants(final List<SageVariant> variants)
    {
        List<VariantData> indels = Lists.newArrayList();
        List<VariantData> candidates = Lists.newArrayList();

        for(SageVariant variant : variants)
        {
            if(!variant.hasLocalPhaseSets())
                continue;

            VariantData adjustedVariant = new VariantData(variant);

            if(variant.isIndel() && adjustedVariant.allowByFilter())
                indels.add(adjustedVariant);

            candidates.add(adjustedVariant);
        }

        if(indels.isEmpty())
            return;

        if(indels.size() > 1)
            Collections.sort(indels, new IndelScoreSorter());

        if(candidates.size() > 1)
            Collections.sort(candidates, new VariantFlankStartSorter());

        while(!indels.isEmpty())
        {
            VariantData indel = indels.get(0);
            indels.remove(0);

            List<VariantData> dedupGroup = findDedupGroup(indel, candidates);

            if(dedupGroup.isEmpty())
                continue;

            // de-dup this group
            dedupIndelGroup(indel, dedupGroup);

            // remove any filtered variants from the main candidate list and the indel list
            for(VariantData variant : dedupGroup)
            {
                candidates.remove(variant);

                if(variant.Variant.isIndel())
                    indels.remove(variant);
            }
        }
    }

    private static final int INDEL_DEDUP_PHASED_DIST_THRESHOLD = 60;
    private static final int LARGE_DEDUP_GROUP_SIZE = 6;
    private static final int LARGE_DEDUP_SELECT_MAX = 3;
    private static final int INDEL_DEDUP_MAX_ITERATIONS = 100;

    private List<VariantData> findDedupGroup(final VariantData indel, final List<VariantData> candidates)
    {
        List<VariantData> dedupGroup = Lists.newArrayList();

        // look for phased variants within the INDELs flanks or with a low maxEdgeDistance
        for(VariantData variant : candidates)
        {
            if(variant == indel)
                continue;

            if(variant.positionEnd() < indel.FlankPosStart - INDEL_DEDUP_PHASED_DIST_THRESHOLD)
                continue;

            if(variant.position() > indel.FlankPosEnd + INDEL_DEDUP_PHASED_DIST_THRESHOLD)
                break;

            if(!indel.Variant.hasMatchingLps(variant.Variant.localPhaseSets()))
                continue;

            if(isDedupCandidate(indel, variant, false))
            {
                dedupGroup.add(variant);
            }
        }

        return dedupGroup;
    }

    private void dedupIndelGroup(final VariantData indel, final List<VariantData> dedupGroup)
    {
        SG_LOGGER.trace("indel({}) with {} other variants", indel, dedupGroup.size());

        // skip if a non-INDEL scores higher than the top indel
        if(dedupGroup.stream().filter(x -> !x.Variant.isIndel()).anyMatch(x -> x.indelScore() > indel.IndelScore))
            return;

        int indelPosition = indel.position();

        List<VariantData> overlappedIndels = Lists.newArrayList();

        // any DEL overlapping the main INDEL can be de-duped immediately
        boolean hasPassingVariant = false;

        int index = 0;
        while(index < dedupGroup.size())
        {
            VariantData variant = dedupGroup.get(index);

            hasPassingVariant |= variant.allowByFilter();

            if(variant.Variant.isDelete() && positionsOverlap(indelPosition, indelPosition, variant.position(), variant.positionEnd()))
            {
                dedupGroup.remove(index);
                overlappedIndels.add(variant);
            }
            else
            {
                ++index;
            }
        }

        // must be at least one passing variant in the group
        if(!hasPassingVariant)
            return;

        String refBases = mRefGenome.getBaseString(indel.Variant.chromosome(), indel.FlankPosStart, indel.FlankPosEnd);

        if(refBases == null || refBases.isEmpty())
            return;

        String indelCoreFlankBases = indel.ReadCounter.readContext().readBases();

        mGroupIterations = 0;

        List<VariantData> dedupedVariants = Lists.newArrayListWithCapacity(dedupGroup.size());

        if(!checkDedupCombinations(indel, dedupGroup, dedupedVariants, indelCoreFlankBases, refBases, indel.FlankPosStart, indel.FlankPosEnd))
        {
            overlappedIndels.forEach(x -> markAsDedup(x.Variant));
            dedupGroup.addAll(overlappedIndels); // add back in so they're removed from further consideration
            return;
        }

        dedupGroup.addAll(overlappedIndels);
        dedupedVariants.addAll(overlappedIndels);
        dedupGroup.add(indel); // so it can be rescued if required

        List<VariantData> nonDedupedVariants = dedupGroup.stream().filter(x -> !dedupedVariants.contains(x)).collect(Collectors.toList());

        for(VariantData variant : dedupGroup)
        {
            if(dedupedVariants.contains(variant))
            {
                // only de-dup variants which fall within the INDEL's bounds
                if(!isDedupCandidate(indel, variant, true))
                    continue;

                markAsDedup(variant.Variant);
            }
            else if(recoverFilteredVariant(variant.Variant, nonDedupedVariants))
            {
                variant.Variant.filters().clear();
            }
        }
    }

    private static boolean recoverFilteredVariant(final SageVariant variant, final List<VariantData> nonDedupedVariants)
    {
        if(variant.isPassing())
            return false;

        if(variant.isIndel() && variant.localPhaseSets().size() == 1 && nonDedupedVariants.size() > 1)
        {
            // check LPS conditions for other variants required to recover this variant
            double maxMatchedLpsReadCountPercent = 0;
            int variantLps = variant.localPhaseSets().get(0);

            for(VariantData otherVariant : nonDedupedVariants)
            {
                if(otherVariant.Variant == variant)
                    continue;

                int otherVarLpsReadCountTotal = 0;
                int otherVarMatchedLpsReadCount = 0;

                for(int i = 0; i < otherVariant.Variant.localPhaseSets().size(); ++i)
                {
                    int lpsReadCount = otherVariant.Variant.localPhaseSetCounts().get(i);
                    otherVarLpsReadCountTotal += lpsReadCount;

                    if(otherVariant.Variant.localPhaseSets().get(i) == variantLps)
                    {
                        otherVarMatchedLpsReadCount = lpsReadCount;
                    }
                }

                double matchedLpsReadCountPercent = otherVarMatchedLpsReadCount / (double)otherVarLpsReadCountTotal;

                maxMatchedLpsReadCountPercent = max(maxMatchedLpsReadCountPercent, matchedLpsReadCountPercent);
            }

            if(maxMatchedLpsReadCountPercent < INDEL_DEDUP_MIN_MATCHED_LPS_PERCENT)
                return false;
        }

        return variant.filters().stream().noneMatch(x -> SoftFilter.GERMLINE_FILTERS.contains(x));
    }

    private void markAsDedup(final SageVariant variant)
    {
        // mark only otherwise passing
        variant.filters().add(DEDUP_INDEL);
    }

    private boolean isDedupCandidate(final VariantData indel, final VariantData variant, boolean requireCoreEndInclusion)
    {
        if(requireCoreEndInclusion)
        {
            if(positionsWithin(variant.position(), variant.CorePosEnd, indel.FlankPosStart, indel.FlankPosEnd))
                return true;
        }
        else
        {
            if(positionWithin(variant.position(), indel.FlankPosStart, indel.FlankPosEnd))
                return true;
        }

        if(variant.ReadCounter.readEdgeDistance().maxAltDistanceFromEdge() < mReadEdgeDistanceThreshold)
            return true;

        return false;
    }

    private boolean checkDedupCombinations(
            final VariantData indel, final List<VariantData> dedupGroup, final List<VariantData> dedupedVariants,
            final String indelCoreFlankBases, final String refBases, final int refPosStart, final int refPosEnd)
    {
        // check for combinations of variants which together explain the sum of differences
        // variants are taken from within the flank positions of the indel

        // first test just adding the indel back to the ref and checking for a match
        if(checkDedupCombination(
                indel, dedupGroup, Collections.emptyList(), dedupedVariants, indelCoreFlankBases, refBases, refPosStart, refPosEnd))
        {
            return true;
        }

        // then check all variants in increasing sized combinations, always plus
        Set<VariantData> variantSet = Sets.newHashSet(dedupGroup);

        int maxSelectionCount = dedupGroup.size() <= LARGE_DEDUP_GROUP_SIZE ? dedupGroup.size() : LARGE_DEDUP_SELECT_MAX;

        for(int i = 1; i <= maxSelectionCount; ++i)
        {
            Set<Set<VariantData>> subsets = Sets.combinations(variantSet, i);

            List<List<VariantData>> variantLists = Lists.newArrayList();

            for(Set<VariantData> subset : subsets)
            {
                variantLists.add(subset.stream().sorted().collect(Collectors.toList()));
            }

            Collections.sort(variantLists, new VariantComboSetSorter());

            for(List<VariantData> variantList : variantLists)
            {
                if(checkDedupCombination(
                        indel, dedupGroup, variantList, dedupedVariants, indelCoreFlankBases, refBases, refPosStart, refPosEnd))
                {
                    return true;
                }
            }

            if(mGroupIterations >= INDEL_DEDUP_MAX_ITERATIONS)
            {
                SG_LOGGER.debug("indel({}) deduped all {} variants, at iteration limit({} select=={})",
                        indel, dedupGroup.size(), mGroupIterations, i);

                // add them all, so only the INDEL will be kept plus any outside the flanks with high enough max edge distance
                dedupedVariants.addAll(dedupGroup);
                return true;
            }
        }

        return false;
    }

    private boolean checkDedupCombination(
            final VariantData indel, final List<VariantData> allVariants, final List<VariantData> selectedVariants, final List<VariantData> dedupedVariants,
            final String indelCoreFlankBases, final String refBases, final int refPosStart, final int refPosEnd)
    {
        ++mGroupIterations;

        List<SimpleVariant> testVariants = Lists.newArrayList(indel.Variant.variant());
        selectedVariants.forEach(x -> testVariants.add(x.Variant.variant()));

        String netAltBases = buildAltBasesString(refBases, refPosStart, refPosEnd, testVariants);

        boolean matched = indelCoreFlankBases.length() >= netAltBases.length() ?
                indelCoreFlankBases.contains(netAltBases) : netAltBases.contains(indelCoreFlankBases);

        if(matched)
        {
            allVariants.stream().filter(x -> !selectedVariants.contains(x)).forEach(x -> dedupedVariants.add(x));
            return true;
        }

        return false;
    }

    @VisibleForTesting
    public static String buildAltBasesString(
            final String refBases, final int refPosStart, final int refPosEnd, final List<SimpleVariant> variants)
    {
        Collections.sort(variants, new VariantReversePositionSorter());

        String altBases = refBases;

        // add variants into the ref bases from right to left so the earlier positions remain unafffected by the added alts
        for(SimpleVariant variant : variants)
        {
            if(!positionWithin(variant.position(), refPosStart, refPosEnd))
                continue;

            int variantRelativePosition = variant.position() - refPosStart;
            String newAltBases = altBases.substring(0, variantRelativePosition);

            newAltBases += variant.alt();

            int refBaseGap = variant.isInsert() ? 1 : variant.ref().length();

            if(variantRelativePosition + refBaseGap <= altBases.length())
                newAltBases += altBases.substring(variantRelativePosition + refBaseGap);

            altBases = newAltBases;
        }

        return altBases;
    }

    public static class VariantFlankStartSorter implements Comparator<VariantData>
    {
        public int compare(final VariantData first, final VariantData second)
        {
            if(first.FlankPosStart == second.FlankPosStart)
                return 0;

            return first.FlankPosStart < second.FlankPosStart ? -1 : 1;
        }
    }

    public static class VariantReversePositionSorter implements Comparator<SimpleVariant>
    {
        public int compare(final SimpleVariant first, final SimpleVariant second)
        {
            if(first.Position == second.Position)
                return 0;

            return first.Position > second.Position ? -1 : 1;
        }
    }

    public static class IndelScoreSorter implements Comparator<VariantData>
    {
        public int compare(final VariantData first, final VariantData second)
        {
            if(first.IndelScore == second.IndelScore)
                return 0;

            return first.IndelScore > second.IndelScore ? -1 : 1;
        }
    }

    public static class VariantComboSetSorter implements Comparator<List<VariantData>>
    {
        public int compare(final List<VariantData> first, final List<VariantData> second)
        {
            for(int i = 0; i < first.size(); ++i)
            {
                int compare = first.get(i).compareTo(second.get(i));

                if(compare != 0)
                    return compare;
            }

            return 0;
        }
    }

    private static final int INDEL_LENGTH_FACTOR = 8;
    private static final int MIN_EVENTS_FACTOR = 12;

    private class VariantData implements Comparable<VariantData>
    {
        public final SageVariant Variant;
        public final ReadContextCounter ReadCounter;
        public final int IndelScore;
        public final int CorePosStart;
        public final int CorePosEnd;
        public final int FlankPosStart;
        public final int FlankPosEnd;

        public VariantData(final SageVariant variant)
        {
            Variant = variant;
            ReadCounter = variant.tumorReadCounters().get(0);

            IndelScore = Variant.isIndel() ? indelScore() : 0;

            // flank positions are estimate since they aren't aware of other variants in their core and flanks
            VariantReadContext readContext = ReadCounter.readContext();

            FlankPosStart = readContext.AlignmentStart;
            FlankPosEnd = readContext.AlignmentEnd;

            CorePosStart = readContext.CorePositionStart;
            CorePosEnd = readContext.CorePositionEnd;
        }

        public String ref()
        {
            return Variant.ref();
        }
        public String alt()
        {
            return Variant.alt();
        }

        public int position() { return Variant.position(); }
        public int positionEnd() { return Variant.variant().positionEnd(); }

        public boolean allowByFilter()
        {
            return Variant.isPassing();

            /* no inclusion of germline filtered variants any longer
            if(Variant.isPassing())
                return true;

            // must only have the germline filters below
            return Variant.filters().stream()
                    .allMatch(x -> x.equals(MAX_GERMLINE_VAF) || x.equals(MAX_GERMLINE_RELATIVE_VAF) || x.equals(MAX_GERMLINE_ALT_SUPPORT));
            */
        }

        public int indelScore()
        {
            return (Variant.isIndel() ? ReadCounter.variant().indelLengthAbs() * INDEL_LENGTH_FACTOR : 0)
                        + ReadCounter.readEdgeDistance().maxAltDistanceFromEdge()
                        - MIN_EVENTS_FACTOR * ReadCounter.minNumberOfEvents();
        }

        @Override
        public int compareTo(final VariantData other)
        {
            if(position() == other.position())
            {
                if(ref().length() == other.ref().length())
                    return 0;

                return ref().length() < other.ref().length() ? -1 : 1;
            }

            return position() < other.position() ? -1 : 1;
        }

        public String toString()
        {
            return format("var(%s:%d %s>%s) corePos(%d - %d) flankPos(%d - %d) distFromEdge(%d) score(%d) filters(%s)",
                Variant.chromosome(), Variant.position(), Variant.ref(), Variant.alt(), CorePosStart, CorePosEnd, FlankPosStart,
                    FlankPosEnd, ReadCounter.readEdgeDistance().maxAltDistanceFromEdge(), IndelScore, Variant.filtersStr());
        }
    }
}
