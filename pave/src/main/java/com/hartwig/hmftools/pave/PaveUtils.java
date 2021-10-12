package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.pave.PaveConstants.GENE_UPSTREAM_DISTANCE;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public final class PaveUtils
{
    public static boolean withinTransRange(final TranscriptData transData, int posStart, int posEnd)
    {
        if(transData.posStrand())
            return positionsOverlap(posStart, posEnd, transData.TransStart - GENE_UPSTREAM_DISTANCE, transData.TransEnd);
        else
            return positionsOverlap(posStart, posEnd, transData.TransStart, transData.TransEnd + GENE_UPSTREAM_DISTANCE);
    }

    public static VariantData createRightAlignedVariant(final VariantData variant, final RefGenomeInterface refGenome)
    {
        if(variant.microhomology().isEmpty() || variant.microhomology().equals(".") || variant.repeatCount() == 0)
            return null;

        int mcLength = variant.microhomology().length();
        int shiftLength = mcLength * (variant.repeatCount() - 1);

        // eg pos 100 AG>A with G repeated 5 times
        int newPosition = variant.Position + shiftLength;

        String newRef = refGenome.getBaseString(variant.Chromosome, newPosition, newPosition + variant.Ref.length() - 1);

        String newAlt = variant.isInsert() ? variant.Alt : newRef.substring(0, 1);

        VariantData raVariant = new VariantData(variant.Chromosome, newPosition, newRef, newAlt);
        raVariant.setRefData(raVariant.refData());
        raVariant.setContext(variant.context());

        // due to complexity, skip trying to phase realigned variants
        raVariant.setVariantDetails(VariantData.NO_LOCAL_PHASE_SET, variant.microhomology(), variant.repeatCount());
        return raVariant;
    }

}
