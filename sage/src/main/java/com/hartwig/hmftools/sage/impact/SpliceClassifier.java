package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_ACCEPTOR_BASES;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_DONOR_EXON_BASES;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_REGION_DISTANCE;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.VariantConsequence;

public class SpliceClassifier
{
    private final RefGenomeInterface mRefGenome;

    public SpliceClassifier(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public static boolean isWithinSpliceRegion(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        int exonStartSpliceRegionStart = exon.Start - SPLICE_REGION_DISTANCE;
        int exonStartSpliceRegionEnd;

        int exonEndSpliceRegionStart;
        int exonEndSpliceRegionEnd = exon.End + SPLICE_REGION_DISTANCE;

        if(transData.Strand == POS_STRAND)
        {
            exonStartSpliceRegionEnd = exon.Start + SPLICE_ACCEPTOR_BASES;
            exonEndSpliceRegionStart = exon.End - SPLICE_DONOR_EXON_BASES;
        }
        else
        {
            exonStartSpliceRegionEnd = exon.Start + SPLICE_DONOR_EXON_BASES;
            exonEndSpliceRegionStart = exon.End - SPLICE_ACCEPTOR_BASES;
        }

        if(exonStartSpliceRegionStart >= variant.Position && variant.Position <= exonStartSpliceRegionEnd)
            return true;
        else if(exonEndSpliceRegionStart >= variant.Position && variant.Position <= exonEndSpliceRegionEnd)
            return true;
        else
            return false;
    }

    public VariantTransImpact classifyVariant(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        // TEMP:
        return new VariantTransImpact(transData, VariantConsequence.SPLICE_REGION_VARIANT);
        /*
                if(transData.Strand == POS_STRAND)
        {
            if(position < exon.Start)
            {
                if(position < exon.Start - SPLICE_REGION_DISTANCE)
                    return false;

                if(position < exon.Start - SPLICE_ACCEPTOR_BASES)
            }
            spliceRegionStart = exon.Start - SPLICE_DONOR_BASES;

        }
        else
        {

        }

         */

        //return null;
    }
}
