package com.hartwig.hmftools.isofox.fusion;

public class ChrGeneCollectionPair
{
    public final String Chromosome;
    public final int GeneCollectionId;

    public ChrGeneCollectionPair(final String chromosome, final int geneCollectionId)
    {
        Chromosome = chromosome;
        GeneCollectionId = geneCollectionId;
    }

    public String toString() { return String.format("%s_%d", Chromosome, GeneCollectionId); }

    public static ChrGeneCollectionPair from(final String chrGenePair)
    {
        final String[] items = chrGenePair.split("_");
        return items.length == 2 ? new ChrGeneCollectionPair(items[0], Integer.parseInt(items[1])) : null;
    }

    public boolean equals(final ChrGeneCollectionPair other)
    {
        return Chromosome.equals(other.Chromosome) && GeneCollectionId == other.GeneCollectionId;
    }
}
