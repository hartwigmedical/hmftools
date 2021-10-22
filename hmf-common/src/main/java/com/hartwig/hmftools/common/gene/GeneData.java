package com.hartwig.hmftools.common.gene;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;

public class GeneData
{
    public final String GeneId; // aka StableId
    public final String GeneName;
    public final String Chromosome;
    public final byte Strand;
    public final int GeneStart;
    public final int GeneEnd;
    public final String KaryotypeBand;

    private String mSynonyms;

    public static final String SYNONYM_DELIM = ";";

    public GeneData(
            final String geneId, final String geneName, final String chromosome, final byte strand,
            final int geneStart, final int geneEnd, final String karyotypeBand)
    {
        GeneId = geneId;
        GeneName = geneName;
        Chromosome = chromosome;
        Strand = strand;
        GeneStart = geneStart;
        GeneEnd = geneEnd;
        KaryotypeBand = karyotypeBand;
        mSynonyms = "";
    }

    public boolean forwardStrand() { return Strand == POS_STRAND; }
    public boolean reverseStrand() { return Strand == NEG_STRAND; }

    public void setSynonyms(final String synonyms) { mSynonyms = synonyms; }
    public boolean hasSynonym(final String name) { return mSynonyms.contains(name); }
    public String getSynonyms() { return mSynonyms; }

    public int length() { return GeneEnd - GeneStart; }

    public String toString()
    {
        return String.format("%s:%s chr(%s) pos(%d-%d) strand(%d)",
                GeneId, GeneName, Chromosome, GeneStart, GeneEnd, Strand);
    }

}
