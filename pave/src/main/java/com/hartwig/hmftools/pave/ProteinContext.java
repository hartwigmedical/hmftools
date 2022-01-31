package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.PaveConstants.DELIM;

import java.util.List;
import java.util.StringJoiner;

import org.apache.commons.compress.utils.Lists;

public class ProteinContext
{
    public String RefCodonBases; // coding bases rounded expanded to cover whole codons
    public String AltCodonBases; // as above but with ref swapped for alt
    public List<int[]> RefCodonsRanges; // the range of the ref codon bases, use for phasing variants

    public String AltCodonBasesComplete; // alt codons plus any downstream refs to make a complete codon, also extended for inframe
    public String RefCodonBasesExtended; // ref codon bases extended for inframe variants until diff AA is found

    public int CodonIndex; // amino acid index of ref codon, corresponds to the coding context CodingBase

    public String RefAminoAcids;
    public String AltAminoAcids;

    // strips off any ref codon present in both
    public boolean ExtraUpstreamCodon;
    public int[] NetCodonIndexRange;
    public String NetRefAminoAcids;
    public String NetAltAminoAcids;
    public boolean IsDuplication;
    public boolean IsPhased;

    public String Hgvs;

    public ProteinContext()
    {
        RefCodonBases = "";
        AltCodonBases = "";
        AltCodonBasesComplete = "";
        RefCodonBasesExtended = "";
        RefCodonsRanges = Lists.newArrayList();

        CodonIndex = 0;
        RefAminoAcids = "";
        AltAminoAcids = "";

        NetCodonIndexRange = new int[] {0, 0};
        NetRefAminoAcids = "";
        NetAltAminoAcids = "";
        IsDuplication = false;
        ExtraUpstreamCodon = false;
        IsPhased = false;
        Hgvs = "";
    }

    public boolean hasCodingBases() { return !RefCodonBases.isEmpty(); }
    public boolean hasAminoAcids() { return !RefAminoAcids.isEmpty(); }

    public boolean hasProteinChange() { return !RefAminoAcids.equals(AltAminoAcids); }

    public boolean validRefCodon() { return !RefCodonBases.isEmpty() && isCodonMultiple(RefCodonBases.length()); }
    public boolean validAltCodon() { return isCodonMultiple(AltCodonBasesComplete.length()); }

    public int refCodingBaseStart() { return refCodingBasePosition(SE_START); }
    public int refCodingBaseEnd() { return refCodingBasePosition(SE_END); }

    public int refCodingBasePosition(int seIndex)
    {
        if(RefCodonsRanges.isEmpty())
            return -1;

        if(seIndex == SE_START)
            return RefCodonsRanges.get(0)[SE_START];
        else
            return RefCodonsRanges.get(RefCodonsRanges.size() - 1)[SE_END];
    }

    public static String csvHeader()
    {
        return "HgvsProtein,RefCodonBases,AltCodonBases,CodonIndex,RefAA,AltAA,NetAARange,NetRefAA,NetAltAA,IsDup";
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
