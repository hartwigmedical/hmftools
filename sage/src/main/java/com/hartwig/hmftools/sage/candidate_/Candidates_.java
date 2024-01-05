package com.hartwig.hmftools.sage.candidate_;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.SimpleVariantComparator;
import com.hartwig.hmftools.sage.select.TierSelector;

import org.apache.commons.compress.utils.Lists;

/**
 * A collection of Candidates that we can add alt contexts to.
 * <p>
 * mSavedCandidates - is really map from AltContext_ as SimpleVariant to Candidate_s.
 */
public class Candidates_
{
    private final List<VariantHotspot> mHotspots;
    private final List<BaseRegion> mPanel;
    private final List<BaseRegion> mHighConfidence;
    private Map<VariantHotspot,List<Candidate_>> mCandidateMap;
    private final SortedMap<SimpleVariant, List<Candidate_>> mSavedCandidates;

    public Candidates_(final List<VariantHotspot> hotspots, final List<BaseRegion> panel, final List<BaseRegion> highConfidence)
    {
        mHotspots = hotspots;
        mPanel = panel;
        mHighConfidence = highConfidence;

        mCandidateMap = null;
        mSavedCandidates = Maps.newTreeMap(new SimpleVariantComparator());
    }

    /**
     * Converts list of alt contexts to candidates and save them.
     */
    public void addSingleSample(final Collection<AltContext_> altContexts)
    {
        final TierSelector tierSelector = new TierSelector(mHotspots, mPanel, mHighConfidence);
        for (AltContext_ altContext : altContexts)
        {
            Candidate_ candidate = Candidate_.fromAltContext(tierSelector.tier(altContext), altContext);
            saveCandidate(candidate);
        }
    }

    /**
     * Creates a Candidate_ from each alt context and add them to mCandidateMap, except when there already exists a candidate with same
     * variant and core.
     * In this latter case we simply swap out the existing candidate with the new alt context if the new alt context has more support.
     */
    public void addOfMultipleSamples(final Collection<AltContext_> altContexts)
    {
        if(mCandidateMap == null)
            mCandidateMap = Maps.newHashMap();

        final TierSelector tierSelector = new TierSelector(mHotspots, mPanel, mHighConfidence);
        final SimpleVariantComparator variantComparator = new SimpleVariantComparator();
        for (AltContext_ altContext : altContexts)
        {
            Candidate_ newCandidate = Candidate_.fromAltContext(tierSelector.tier(altContext), altContext);
            List<Candidate_> candidates = mCandidateMap.get(altContext);
            if (candidates == null)
            {
                candidates = Lists.newArrayList();
                mCandidateMap.put(altContext, candidates);
            }

            // Does there already exist a candidate that has the same variant and core string?
            Candidate_ matchingCandidate = candidates.stream()
                    .filter(candidate -> variantComparator.compare(candidate.variant(), newCandidate.variant()) == 0)
                    .filter(candidate -> candidate.readContext().coreString().equals(newCandidate.readContext().coreString()))
                    .findFirst().orElse(null);

            if (matchingCandidate != null)
            {
                // Swap them out if new altContext has more support.
                matchingCandidate.update(altContext);
            }
            else
            {
                candidates.add(newCandidate);
            }
        }
    }

    /**
     * Returns the saved candidates (sorted via SimpleVariantComparator), with mCandidateMap values absorbed in.
     * @param restrictedPositions Filters the saved candidates to ones whose positions are contained in this.
     * @return Filtered and sorted saved candidates.
     */
    public List<Candidate_> candidates(final Set<Integer> restrictedPositions)
    {
        if (mCandidateMap != null)
            mCandidateMap.values().stream().flatMap(List::stream).forEach(this::saveCandidate);

        if (restrictedPositions != null && !restrictedPositions.isEmpty())
        {
            Stream<Candidate_> allCandidates = mSavedCandidates.values().stream().flatMap(List::stream);
            List<Candidate_> restrictedCandidates = allCandidates.filter(candidate -> restrictedPositions.contains(candidate.position())).collect(Collectors.toList());

            mSavedCandidates.clear();
            restrictedCandidates.forEach(this::saveCandidate);
        }

        return mSavedCandidates.values().stream().flatMap(List::stream).collect(Collectors.toList());
    }

    /**
     * Add candidate to saved candidates.
     */
    private void saveCandidate(final Candidate_ candidate)
    {
        SimpleVariant variantKey = candidate.variant();
        List<Candidate_> candidates = mSavedCandidates.get(variantKey);
        if (candidates == null)
        {
            candidates = Lists.newArrayList();
            mSavedCandidates.put(variantKey, candidates);
        }

        candidates.add(candidate);
    }
}
