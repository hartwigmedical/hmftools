package com.hartwig.hmftools.common.gene;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_DELIM;

import java.util.StringJoiner;

public class TranscriptAminoAcids
{
    public final String GeneId;
    public final String GeneName;
    public final String TransName;
    public final String AminoAcids;

    public TranscriptAminoAcids(final String geneId, final String geneName, final String transName, final String aminoAcids)
    {
        GeneId = geneId;
        GeneName = geneName;
        TransName = transName;
        AminoAcids = aminoAcids;
    }

    public static String csvHeader() { return "GeneId,GeneName,TransName,AminoAcids"; }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(ENSEMBL_DELIM);
        sj.add(GeneId);
        sj.add(GeneName);
        sj.add(TransName);
        sj.add(AminoAcids);
        return sj.toString();
    }

    public int baseLength() { return AminoAcids.length() + 1; }

}
