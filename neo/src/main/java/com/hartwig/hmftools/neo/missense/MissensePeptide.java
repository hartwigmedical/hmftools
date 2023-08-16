package com.hartwig.hmftools.neo.missense;

import static java.lang.String.format;

public class MissensePeptide
{
    public String GeneId;
    public String GeneName;
    public String TransName;

    public int Position; // of the missense variant
    public int CodonIndex; // of the codon affected by the variant
    public String Context;
    public char RefBase;
    public char AltBase;

    public String Peptide;
    public String UpFlank;
    public String DownFlank;
    public int PeptideStartIndex; // codon index for start of peptide

    public MissensePeptide(final String geneId, final String geneName, final String transName)
    {
        GeneId = geneId;
        GeneName = geneName;
        TransName = transName;
        Peptide = "";
        UpFlank = "";
        DownFlank = "";
        PeptideStartIndex = -1;
    }

    public String toString()
    {
        return format("%s var(%d %c>%c) codon(%d start=%d) peptide(%s) flanks(up=%s down=%s)",
                GeneName, Position, RefBase, AltBase, CodonIndex, PeptideStartIndex, Peptide, UpFlank, DownFlank);
    }
}
