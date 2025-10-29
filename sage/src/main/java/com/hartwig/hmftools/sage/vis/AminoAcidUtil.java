package com.hartwig.hmftools.sage.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.sage.vis.AminoAcidVariant.AminoAcidSubstitution;
import com.hartwig.hmftools.sage.vis.GeneRegionViewModel.AminoAcidViewModel_;
import com.hartwig.hmftools.sage.vis.GeneRegionViewModel.IntronicRegionViewModel;

public final class AminoAcidUtil
{
    private AminoAcidUtil() {}

    public static List<BaseRegion> getCodingRegions(final TranscriptData transcriptExons)
    {
        List<BaseRegion> codingRegions = Lists.newArrayList();
        BaseRegion transcriptCodingRegion = new BaseRegion(transcriptExons.CodingStart, transcriptExons.CodingEnd);
        List<ExonData> exons = transcriptExons.exons();
        for(ExonData exon : exons)
        {
            BaseRegion exonRegion = new BaseRegion(exon.Start, exon.End);
            if(!transcriptCodingRegion.overlaps(exonRegion))
                continue;

            int codingStart = max(transcriptCodingRegion.start(), exonRegion.start());
            int codingEnd = min(transcriptCodingRegion.end(), exonRegion.end());
            BaseRegion codingRegion = new BaseRegion(codingStart, codingEnd);
            codingRegions.add(codingRegion);
        }

        return codingRegions;
    }

    public static List<GeneRegionViewModel> getGeneRegions(
            final TranscriptData transcriptExons, final TranscriptAminoAcids transcriptAminoAcids, final List<AminoAcidVariant> variants)
    {
        List<GeneRegionViewModel> geneRegions = Lists.newArrayList();
        boolean posStrand = transcriptExons.Strand == (byte) 1;
        String aminoAcids = transcriptAminoAcids.AminoAcids;
        List<BaseRegion> codingRegions = getCodingRegions(transcriptExons);
        int nucIdx = 0;
        int aaIdx = posStrand ? 0 : aminoAcids.length() - 1;
        BaseRegion prevCodingRegion = null;
        for(int i = 0; i < codingRegions.size(); i++)
        {
            BaseRegion codingRegion = codingRegions.get(i);
            if(prevCodingRegion != null)
                geneRegions.add(new IntronicRegionViewModel(new BaseRegion(prevCodingRegion.end() + 1, codingRegion.start() - 1)));

            prevCodingRegion = codingRegion;

            int posStart = codingRegion.start();
            while(posStart <= codingRegion.end())
            {
                int posEnd = min(posStart + 2 - nucIdx, codingRegion.end());
                nucIdx = (posEnd - posStart + 1 + nucIdx) % 3;
                char refAcid = aminoAcids.charAt(aaIdx);
                int aaPos = aaIdx + 1;
                BaseRegion aaRegion = new BaseRegion(posStart, posEnd);
                char altAcid = refAcid;
                for(AminoAcidVariant variant : variants)
                {
                    if(variant instanceof AminoAcidSubstitution)
                    {
                        AminoAcidSubstitution missenseVariant = (AminoAcidSubstitution) variant;
                        if(missenseVariant.aminoAcidPos() != aaPos)
                            continue;

                        altAcid = missenseVariant.alt();
                    }
                    else
                    {
                        throw new RuntimeException(format("Cannot handle variant of type: %s", variant.getClass().getSimpleName()));
                    }
                }

                geneRegions.add(new AminoAcidViewModel_(aaRegion, aaPos, refAcid, altAcid));
                if(nucIdx == 0)
                    aaIdx += posStrand ? 1 : -1;

                posStart = posEnd + 1;
            }
        }

        if(posStrand && aaIdx != aminoAcids.length() || !posStrand && aaIdx != -1)
            throw new RuntimeException("Transcript amino acid count does not match coding length.");

        return geneRegions;
    }
}
