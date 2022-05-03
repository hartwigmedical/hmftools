package com.hartwig.hmftools.isofox.fusion;

public class ChimericStats
{
    public int ChimericJunctions;
    public int LocalInterGeneFrags;
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
        LocalInterGeneFrags += other.LocalInterGeneFrags;
        CandidateRealignFrags += other.CandidateRealignFrags;
        Translocations += other.Translocations;
        Inversions += other.Inversions;
    }

    public void clear()
    {
        ChimericJunctions = 0;
        LocalInterGeneFrags = 0;
        CandidateRealignFrags = 0;
        Translocations = 0;
        Inversions = 0;
    }

    public String toString()
    {
        return String.format("junc=%d locInterGene=%d candRealgn=%d bnd=%d",
                ChimericJunctions, LocalInterGeneFrags, CandidateRealignFrags, Translocations);
    }

}
