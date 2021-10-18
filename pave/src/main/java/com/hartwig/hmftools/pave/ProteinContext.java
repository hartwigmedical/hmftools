package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.PaveConstants.DELIM;

import java.util.StringJoiner;

public class ProteinContext
{
    public String RefCodonBases; // coding bases rounded expanded to cover whole codons
    public String AltCodonBases; // as above but with ref swapped for alt
    public String AltCodonBasesComplete; // alt codons plus any subsequent downstream refs to make complete codon(s)
    public int[] RefCodonsRange; // the range of the ref codon bases, use for phasing variants

    public int CodonIndex; // amino acid index of ref codon, corresponds to the coding context CodingBase
    public String RefAminoAcids;
    public String AltAminoAcids;

    // strips off any ref codon present in both
    public int[] NetCodonIndexRange;
    public String NetRefAminoAcids;
    public String NetAltAminoAcids;
    public boolean IsDuplication;

    public String Hgvs;

    public ProteinContext()
    {
        RefCodonBases = "";
        AltCodonBases = "";
        AltCodonBasesComplete = "";
        RefCodonsRange = new int[] {0, 0};

        CodonIndex = 0;
        RefAminoAcids = "";
        AltAminoAcids = "";

        NetCodonIndexRange = new int[] {0, 0};
        NetRefAminoAcids = "";
        NetAltAminoAcids = "";
        IsDuplication = false;
        Hgvs = "";
    }

    public boolean hasCodingBases() { return !RefAminoAcids.isEmpty(); }

    public boolean hasProteinChange() { return !RefAminoAcids.equals(AltAminoAcids); }

    public boolean validRefCodon() { return !RefCodonBases.isEmpty() && isCodonMultiple(RefCodonBases.length()); }
    public boolean validAltCodon() { return isCodonMultiple(AltCodonBasesComplete.length()); }

    public static String csvHeader()
    {
        return "HgvsProtein,RefCodonBases,AltCodonBases,CodonIndex,RefAA,AltAA,NetAARange,NetRefAA,NegAltAA,IsDup";
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
        sj.add(String.format("%d-%d", NetCodonIndexRange[SE_START], NetCodonIndexRange[SE_END]));
        sj.add(NetRefAminoAcids);
        sj.add(NetAltAminoAcids);
        sj.add(String.valueOf(IsDuplication));
        return sj.toString();
    }

    public static String empty() { return ",,,0,,,0-0,,,"; }
}
