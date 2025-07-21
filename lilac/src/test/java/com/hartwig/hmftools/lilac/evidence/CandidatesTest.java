package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_DEPTH_FILTER;
import static com.hartwig.hmftools.lilac.evidence.Candidates.filterCandidates;
import static com.hartwig.hmftools.lilac.evidence.Candidates.filterSequencesByMinSupport;

import static org.junit.Assert.assertEquals;
import static org.mockito.ArgumentMatchers.anyInt;
import static org.mockito.Mockito.doReturn;
import static org.mockito.Mockito.mock;

import java.util.List;
import java.util.NavigableMap;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.hartwig.hmftools.lilac.LilacConstants;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.util.ThrowOnUnstubbed;

import org.junit.Test;

public class CandidatesTest
{
    @Test
    public void testFilterSequencesByMinSupportCandidateFilteredOut()
    {
        LilacConstants.MIN_DEPTH_FILTER = DEFAULT_MIN_DEPTH_FILTER;

        HlaAllele allele = HlaAllele.fromString("A*01:01");
        List<HlaSequenceLoci> initialCandidates = Lists.newArrayList(
                new HlaSequenceLoci(allele, Lists.newArrayList("M", "M")),
                new HlaSequenceLoci(allele, Lists.newArrayList("M", "V")));

        NavigableMap<Integer, Multiset<String>> pooledSeqCountsByLoci = Maps.newTreeMap();
        pooledSeqCountsByLoci.put(0, HashMultiset.create());
        pooledSeqCountsByLoci.put(1, HashMultiset.create());
        pooledSeqCountsByLoci.get(0).setCount("M", 100);
        pooledSeqCountsByLoci.get(1).setCount("M", 100);

        SequenceCount pooledCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(pooledSeqCountsByLoci).when(pooledCounts).seqCountsByLoci();
        doReturn(List.of("M")).when(pooledCounts).getMinEvidenceSequences(anyInt());

        List<Integer> aminoAcidBoundaries = Lists.newArrayList();

        Multiset<String> localAcids = HashMultiset.create();
        localAcids.setCount("M", DEFAULT_MIN_DEPTH_FILTER);
        SequenceCount localCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(localAcids).when(localCounts).get(anyInt());

        List<HlaSequenceLoci> candidates = filterSequencesByMinSupport(initialCandidates, pooledCounts, aminoAcidBoundaries, localCounts);

        assertEquals(1, candidates.size());
        assertEquals("MM", candidates.get(0).sequence());
    }

    @Test
    public void testFilterSequencesByMinSupportCandidateSavedDueToLowDepth()
    {
        LilacConstants.MIN_DEPTH_FILTER = DEFAULT_MIN_DEPTH_FILTER;

        HlaAllele allele = HlaAllele.fromString("A*01:01");
        List<HlaSequenceLoci> initialCandidates = Lists.newArrayList(
                new HlaSequenceLoci(allele, Lists.newArrayList("M", "M")),
                new HlaSequenceLoci(allele, Lists.newArrayList("M", "V")));

        NavigableMap<Integer, Multiset<String>> pooledSeqCountsByLoci = Maps.newTreeMap();
        pooledSeqCountsByLoci.put(0, HashMultiset.create());
        pooledSeqCountsByLoci.put(1, HashMultiset.create());
        pooledSeqCountsByLoci.get(0).setCount("M", 100);
        pooledSeqCountsByLoci.get(1).setCount("M", 100);

        SequenceCount pooledCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(pooledSeqCountsByLoci).when(pooledCounts).seqCountsByLoci();
        doReturn(List.of("M")).when(pooledCounts).getMinEvidenceSequences(anyInt());

        List<Integer> aminoAcidBoundaries = Lists.newArrayList();

        Multiset<String> localHighAcids = HashMultiset.create();
        localHighAcids.setCount("M", DEFAULT_MIN_DEPTH_FILTER);
        Multiset<String> localLowAcids = HashMultiset.create();
        localLowAcids.setCount("M", DEFAULT_MIN_DEPTH_FILTER - 1);
        SequenceCount localCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(localHighAcids).when(localCounts).get(0);
        doReturn(localLowAcids).when(localCounts).get(1);

        List<HlaSequenceLoci> candidates = filterSequencesByMinSupport(initialCandidates, pooledCounts, aminoAcidBoundaries, localCounts);

        assertEquals(2, candidates.size());
    }

    @Test
    public void testFilterCandidatesCandidateFilteredOut()
    {
        LilacConstants.MIN_DEPTH_FILTER = DEFAULT_MIN_DEPTH_FILTER;

        HlaAllele allele = HlaAllele.fromString("A*01:01");
        List<HlaSequenceLoci> initialCandidates = Lists.newArrayList(
                new HlaSequenceLoci(allele, Lists.newArrayList("M", "L", "V")),
                new HlaSequenceLoci(allele, Lists.newArrayList("M", "L", "L")));

        PhasedEvidence phasedEvidence = mock(PhasedEvidence.class, new ThrowOnUnstubbed());
        doReturn(List.of(0, 2)).when(phasedEvidence).getAminoAcidLoci();
        doReturn(ImmutableMap.of("MV", 0)).when(phasedEvidence).getEvidence();
        List<PhasedEvidence> evidence = Lists.newArrayList(phasedEvidence);

        Multiset<String> localAcids = HashMultiset.create();
        localAcids.setCount("M", DEFAULT_MIN_DEPTH_FILTER);
        SequenceCount localCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(localAcids).when(localCounts).get(anyInt());

        List<HlaSequenceLoci> candidates = filterCandidates(initialCandidates, evidence, localCounts);

        assertEquals(1, candidates.size());
        assertEquals("MLV", candidates.get(0).sequence());
    }

    @Test
    public void testFilterCandidatesCandidateSavedDueToLowDepth()
    {
        LilacConstants.MIN_DEPTH_FILTER = DEFAULT_MIN_DEPTH_FILTER;

        HlaAllele allele = HlaAllele.fromString("A*01:01");
        List<HlaSequenceLoci> initialCandidates = Lists.newArrayList(
                new HlaSequenceLoci(allele, Lists.newArrayList("M", "L", "V")),
                new HlaSequenceLoci(allele, Lists.newArrayList("M", "L", "L")));

        PhasedEvidence phasedEvidence = mock(PhasedEvidence.class, new ThrowOnUnstubbed());
        doReturn(List.of(0, 2)).when(phasedEvidence).getAminoAcidLoci();
        doReturn(ImmutableMap.of("MV", 0)).when(phasedEvidence).getEvidence();
        List<PhasedEvidence> evidence = Lists.newArrayList(phasedEvidence);

        Multiset<String> localHighAcids = HashMultiset.create();
        localHighAcids.setCount("M", DEFAULT_MIN_DEPTH_FILTER);
        Multiset<String> localLowAcids = HashMultiset.create();
        localLowAcids.setCount("M", DEFAULT_MIN_DEPTH_FILTER - 1);
        SequenceCount localCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(localHighAcids).when(localCounts).get(0);
        doReturn(localLowAcids).when(localCounts).get(2);

        List<HlaSequenceLoci> candidates = filterCandidates(initialCandidates, evidence, localCounts);

        assertEquals(2, candidates.size());
    }
}
