package com.hartwig.hmftools.geneutils.panelfinder;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

public class GeneExonData
{
    public final String Gene;
    public final int ExonRank;
    public final int ExonStart;
    public final int ExonEnd;

    public GeneExonData(final String gene, final int exonRank, final int exonStart, final int exonEnd)
    {
        Gene = gene;
        ExonRank = exonRank;
        ExonStart = exonStart;
        ExonEnd = exonEnd;
    }

    public String toString() { return format("%s exon(%d: %d-%d)", Gene, ExonRank, ExonStart, ExonEnd); }

    public static String toString(final List<GeneExonData> geneExons)
    {
        return geneExons.stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
    }

}
