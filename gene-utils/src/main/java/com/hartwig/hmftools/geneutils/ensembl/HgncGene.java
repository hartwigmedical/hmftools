package com.hartwig.hmftools.geneutils.ensembl;

public class HgncGene
{
    public final String HgncId;
    public final String EnsemblGeneId;
    public final String Symbol;

    public HgncGene(final String hgncId, final String ensemblGeneId, final String symbol)
    {
        HgncId = hgncId;
        EnsemblGeneId = ensemblGeneId;
        Symbol = symbol;
    }
}
