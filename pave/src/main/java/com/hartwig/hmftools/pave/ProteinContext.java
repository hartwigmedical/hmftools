package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.pave.PaveConstants.DELIM;

import java.util.StringJoiner;

public class ProteinContext
{
    public String RefCodonBases; // coding bases rounded expanded to cover whole codons
    public String AltCodonBases; // as above but with ref swapped for alt

    public int CodonIndex; // first amino acid affected
    public String RefAminoAcids;
    public String AltAminoAcids;

    public ProteinContext()
    {
        RefCodonBases = "";
        AltCodonBases = "";

        CodonIndex = 0;
        RefAminoAcids = "";
        AltAminoAcids = "";
    }

    public boolean hasCodingBases() { return !RefAminoAcids.isEmpty(); }

    public boolean hasProteinChange() { return !RefAminoAcids.equals(AltAminoAcids); }

    public boolean validRefCodon() { return (RefCodonBases.length() % 3) == 0; }

    public String hgvsStr() { return "tbc"; }

    public static String csvHeader()
    {
        return "HgvsProtein,RefCodonBases,AltCodonBases,CodonIndex,RefAA,AltAA";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIM);

        sj.add(hgvsStr());
        sj.add(RefCodonBases);
        sj.add(AltCodonBases);
        sj.add(String.valueOf(CodonIndex));
        sj.add(RefAminoAcids);
        sj.add(AltAminoAcids);
        return sj.toString();
    }

    public static String empty() { return ",,,0,,"; }
}
