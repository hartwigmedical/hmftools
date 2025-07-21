package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_DEPTH_FILTER;

import static org.junit.Assert.assertEquals;
import static org.mockito.ArgumentMatchers.anyInt;
import static org.mockito.ArgumentMatchers.anySet;
import static org.mockito.Mockito.doAnswer;
import static org.mockito.Mockito.doReturn;
import static org.mockito.Mockito.mock;

import java.util.List;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Set;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConstants;
import com.hartwig.hmftools.lilac.evidence.AminoAcid;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.util.ThrowOnUnstubbed;

import org.junit.Test;

public class AminoAcidTest
{
    @Test
    public void testApplyQualFilterSingleAminoAcidFilteredOut()
    {
        LilacConstants.MIN_DEPTH_FILTER = DEFAULT_MIN_DEPTH_FILTER;
        final Set<Integer> actualFilteredIntersect = Sets.newTreeSet();

        List<AminoAcid> fragmentAcids = Lists.newArrayList(
                new AminoAcid(0, "M"),
                new AminoAcid(1, "V"));
        final NavigableMap<Integer, AminoAcid> fragmentAcidsByLoci = Maps.newTreeMap();
        fragmentAcids.forEach(x -> fragmentAcidsByLoci.put(x.locus(), x));

        Fragment fragment = mock(Fragment.class, new ThrowOnUnstubbed());
        doReturn(fragmentAcidsByLoci).when(fragment).aminoAcidsByLoci();
        fragmentAcids.forEach(x -> doReturn(x.acid()).when(fragment).aminoAcid(x.locus()));
        doAnswer(invocation ->
        {
            Set<Integer> filteredLocus = (Set<Integer>) invocation.getArguments()[0];
            actualFilteredIntersect.addAll(filteredLocus);
            return null;
        }).when(fragment).filterAminoAcidsOnLoci(anySet());

        List<String> minEvidenceSequences = Lists.newArrayList("M");
        SequenceCount pooledCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(minEvidenceSequences).when(pooledCounts).getMinEvidenceSequences(anyInt());

        Multiset<String> localAcids = HashMultiset.create();
        localAcids.setCount("M", DEFAULT_MIN_DEPTH_FILTER);
        SequenceCount localCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(localAcids).when(localCounts).get(anyInt());

        AminoAcidQualEnrichment.applyQualFilter(fragment, pooledCounts, localCounts);
        NavigableSet<Integer> expectedFilteredInserect = Sets.newTreeSet(List.of(0));

        assertEquals(expectedFilteredInserect, actualFilteredIntersect);
    }

    @Test
    public void testApplyQualFilterSingleAminoAcidSavedDueToLowDepth()
    {
        LilacConstants.MIN_DEPTH_FILTER = DEFAULT_MIN_DEPTH_FILTER;
        final Set<Integer> actualFilteredIntersect = Sets.newTreeSet();

        List<AminoAcid> fragmentAcids = Lists.newArrayList(
                new AminoAcid(0, "M"),
                new AminoAcid(1, "V"));
        final NavigableMap<Integer, AminoAcid> fragmentAcidsByLoci = Maps.newTreeMap();
        fragmentAcids.forEach(x -> fragmentAcidsByLoci.put(x.locus(), x));

        Fragment fragment = mock(Fragment.class, new ThrowOnUnstubbed());
        doReturn(fragmentAcidsByLoci).when(fragment).aminoAcidsByLoci();
        fragmentAcids.forEach(x -> doReturn(x.acid()).when(fragment).aminoAcid(x.locus()));
        doAnswer(invocation ->
        {
            Set<Integer> filteredLocus = (Set<Integer>) invocation.getArguments()[0];
            actualFilteredIntersect.addAll(filteredLocus);
            return null;
        }).when(fragment).filterAminoAcidsOnLoci(anySet());

        List<String> minEvidenceSequences = Lists.newArrayList("M");
        SequenceCount pooledCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(minEvidenceSequences).when(pooledCounts).getMinEvidenceSequences(anyInt());

        Multiset<String> localHighAcids = HashMultiset.create();
        localHighAcids.setCount("M", DEFAULT_MIN_DEPTH_FILTER);
        Multiset<String> localLowAcids = HashMultiset.create();
        localLowAcids.setCount("M", DEFAULT_MIN_DEPTH_FILTER - 1);
        SequenceCount localCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(localHighAcids).when(localCounts).get(0);
        doReturn(localLowAcids).when(localCounts).get(1);

        AminoAcidQualEnrichment.applyQualFilter(fragment, pooledCounts, localCounts);
        NavigableSet<Integer> expectedFilteredInserect = Sets.newTreeSet(List.of(0, 1));

        assertEquals(expectedFilteredInserect, actualFilteredIntersect);
    }
}
