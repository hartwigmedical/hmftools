package com.hartwig.hmftools.isofox.fusion;

import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

public class ChimericStats
{
    public int ChimericJunctions;
    public int SupplementaryFrags;
    public int LocalInterGeneFrags;
    public int LocalPreGeneFrags;
    public int CandidateRealignFrags;

    public int Translocations;
    public int Inversions;

    public ChimericStats()
    {
        clear();
    }

    public void merge(final ChimericStats other)
    {
        ChimericJunctions += other.ChimericJunctions;
        SupplementaryFrags += other.SupplementaryFrags;
        LocalInterGeneFrags += other.LocalInterGeneFrags;
        CandidateRealignFrags += other.CandidateRealignFrags;
        Translocations += other.Translocations;
        Inversions += other.Inversions;
        LocalPreGeneFrags += other.LocalPreGeneFrags;
    }

    public void clear()
    {
        ChimericJunctions = 0;
        SupplementaryFrags = 0;
        LocalInterGeneFrags = 0;
        CandidateRealignFrags = 0;
        Translocations = 0;
        Inversions = 0;
        LocalPreGeneFrags = 0;
    }

    public String toString()
    {
        return String.format("junc=%d supp=%d locInterGene=%d candRealgn=%d bnd=%d inv=%d preGene=%d",
                ChimericJunctions, SupplementaryFrags, LocalInterGeneFrags, CandidateRealignFrags,
                Translocations, Inversions, LocalPreGeneFrags);
    }

}
