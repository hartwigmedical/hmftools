package com.hartwig.hmftools.common.ensemblcache;

public class EnsemblGeneData
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

    public EnsemblGeneData(
            String geneId, String geneName, String chromosome, byte strand, int geneStart, int geneEnd, String karyotypeBand)
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

    public boolean forwardStrand() { return Strand == 1; }
    public boolean reverseStrand() { return Strand == -1; }

    public void addSynonyms(final String synonyms) { mSynonyms = synonyms; }
    public boolean hasSynonym(final String name) { return mSynonyms.contains(name); }
    public String getSynonyms() { return mSynonyms; }

    public int length() { return GeneEnd - GeneStart; }

    public String toString()
    {
        return String.format("%s:%s chr(%s) pos(%d-%d) strand(%d)",
                GeneId, GeneName, Chromosome, GeneStart, GeneEnd, Strand);
    }

}
