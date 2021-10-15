package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.pave.PaveConstants.DELIM;

import java.util.StringJoiner;

public class ProteinContext
{
    public String RefCodonBases; // coding bases rounded expanded to cover whole codons
    public String AltCodonBases; // as above but with ref swapped for alt
    public int[] RefCodonsRange; // as above but with ref swapped for alt

    public int CodonIndex; // first amino acid affected, corresponds to the coding context CodingBase
    public String RefAminoAcids;
    public String AltAminoAcids;

    public String Hgvs;

    public ProteinContext()
    {
        RefCodonBases = "";
        AltCodonBases = "";
        RefCodonsRange = new int[] {0, 0};

        CodonIndex = 0;
        RefAminoAcids = "";
        AltAminoAcids = "";
        Hgvs = "";
    }

    public boolean hasCodingBases() { return !RefAminoAcids.isEmpty(); }

    public boolean hasProteinChange() { return !RefAminoAcids.equals(AltAminoAcids); }

    public boolean validRefCodon() { return !RefCodonBases.isEmpty() && isCodonMultiple(RefCodonBases.length()); }
    public boolean validAltCodon() { return isCodonMultiple(AltCodonBases.length()); }

    public static String csvHeader()
    {
        return "HgvsProtein,RefCodonBases,AltCodonBases,CodonIndex,RefAA,AltAA";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIM);

        sj.add(Hgvs);
        sj.add(RefCodonBases);
        sj.add(AltCodonBases);
        sj.add(String.valueOf(CodonIndex));
        sj.add(RefAminoAcids);
        sj.add(AltAminoAcids);
        return sj.toString();
    }

    public static String empty() { return ",,,0,,"; }
}
