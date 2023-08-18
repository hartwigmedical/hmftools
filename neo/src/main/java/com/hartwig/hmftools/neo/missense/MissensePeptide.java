package com.hartwig.hmftools.neo.missense;

import static java.lang.String.format;

public class MissensePeptide
{
    public String GeneId;
    public String GeneName;
    public String TransName;

    public int Position; // of the missense variant
    public int CodonIndex;
    public String Context;
    public char RefBase;
    public char AltBase;

    public String Peptide;
    public String UpFlank;
    public String DownFlank;

    public MissensePeptide(final String geneId, final String geneName, final String transName)
    {
        GeneId = geneId;
        GeneName = geneName;
        TransName = transName;
    }

    public String toString()
    {
        return format("%s var(%d %c>%c) codon(%d) peptide(%s) flanks(up=%s down=%s)",
                GeneName, Position, RefBase, AltBase, CodonIndex, Peptide, UpFlank, DownFlank);
    }
}
