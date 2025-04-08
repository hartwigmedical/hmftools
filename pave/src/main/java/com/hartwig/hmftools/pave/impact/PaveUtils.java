package com.hartwig.hmftools.pave.impact;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.pave.PaveConstants.GENE_UPSTREAM_DISTANCE;

import java.util.List;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pave.GeneCacheIndexing;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.VariantData;

import org.jetbrains.annotations.Nullable;

public final class PaveUtils
{
    public static void findVariantImpacts(
            final VariantData variant, final ImpactClassifier impactClassifier, final GeneDataCache geneDataCache)
    {
        findVariantImpacts(variant, impactClassifier, geneDataCache, null);
    }

    public static void findVariantImpacts(
            final VariantData variant, final ImpactClassifier impactClassifier, final GeneDataCache geneDataCache,
            @Nullable final GeneCacheIndexing cacheIndexing)
    {
        boolean processed = false;

        List<GeneData> geneCandidates = geneDataCache.findGenes(variant.Chromosome, variant.Position, variant.EndPosition, cacheIndexing);

        if(!geneCandidates.isEmpty())
        {
            // analyse against each of the genes and their transcripts
            for(GeneData geneData : geneCandidates)
            {
                List<TranscriptData> transDataList =
                        geneDataCache.findTranscripts(geneData.GeneId, variant.Position, variant.EndPosition);

                // non-coding transcripts are skipped for now
                if(transDataList.isEmpty())
                    continue;

                for(TranscriptData transData : transDataList)
                {
                    VariantTransImpact transImpact = impactClassifier.classifyVariant(variant, transData);
                    processed = true;

                    // check right-alignment if the variant has microhomology
                    if(variant.realignedVariant() != null)
                    {
                        VariantTransImpact raTransImpact = impactClassifier.classifyVariant(variant.realignedVariant(), transData);

                        if(raTransImpact != null)
                        {
                            variant.realignedVariant().addImpact(geneData.GeneName, raTransImpact);
                            transImpact = ImpactClassifier.selectAlignedImpacts(transImpact, raTransImpact);
                        }
                    }

                    if(transImpact != null)
                        variant.addImpact(geneData.GeneName, transImpact);
                }
            }
        }

        // ensure all phased variants are cached
        if(!processed && variant.hasLocalPhaseSet())
            impactClassifier.phasedVariants().checkAddVariant(variant);
    }

    public static boolean withinTransRange(final TranscriptData transData, int posStart, int posEnd)
    {
        if(transData.posStrand())
            return positionsOverlap(posStart, posEnd, transData.TransStart - GENE_UPSTREAM_DISTANCE, transData.TransEnd);
        else
            return positionsOverlap(posStart, posEnd, transData.TransStart, transData.TransEnd + GENE_UPSTREAM_DISTANCE);
    }

    public static VariantData createRightAlignedVariant(final VariantData variant, final RefGenomeInterface refGenome)
    {
        if(variant.isBaseChange()) // to be confirmed
            return null;

        if(variant.microhomology().isEmpty() || variant.microhomology().equals("."))
            return null;

        // cannot realign if the ref and alt bases differ, eg T>AA
        if(variant.Ref.charAt(0) != variant.Alt.charAt(0))
            return null;

        // repeat count can only be used where the alt bases match the microhomology and repeat sequence
        // otherwise shift by the microhomology

        /*
        if the microhomology == ins/del sequence and the MH==repeatSeq, then you can extend to the end of the repeat (-1 repeat sequences for the del case)

        instead:
        if the microhomology == ins/del sequence   and the MH==N*repeatSeq, then you can extend to the end of the repeat ( -N repeat sequences for the del case)
        where N = any multiple of the repeat sequence (in this example N=2)
        */

        String altBases = variant.isDeletion() ? variant.Ref.substring(1) : variant.Alt.substring(1);

        int mcLength = variant.microhomology().length();

        int shiftCount;

        if(variant.microhomology().equals(altBases) && variant.repeatCount() > 0)
        {
            int repeatCount = 0;

            if(variant.repeatSequence().equals(altBases))
            {
                repeatCount = variant.repeatCount();
            }
            else
            {
                // test whether MH = N * repeatSequence
                int repeatSeqLen = variant.repeatSequence().length();

                int nCount = mcLength / repeatSeqLen;

                String testMH = "";

                for(int i = 0; i < nCount; ++i)
                {
                    testMH += variant.repeatSequence();
                }

                if(testMH.equals(variant.microhomology()))
                {
                    repeatCount = variant.repeatCount() / nCount;
                }
            }

            if(repeatCount >= 1)
                shiftCount = variant.isInsert() ? repeatCount : repeatCount - 1;
            else
                shiftCount = 1;
        }
        else
        {
            shiftCount = 1;
        }

        int shiftLength = mcLength * shiftCount;

        // eg pos 100 AG>A with G repeated 5 times
        int newPosition = variant.Position + shiftLength;

        String newRef = refGenome.getBaseString(variant.Chromosome, newPosition, newPosition + variant.Ref.length() - 1);

        String newAlt = variant.isInsert() ?
                newRef.substring(0, 1) + variant.Alt.substring(1) : newRef.substring(0, 1);

        VariantData raVariant = new VariantData(variant.Chromosome, newPosition, newRef, newAlt);
        raVariant.setContext(variant.context());

        // due to complexity, skip trying to phase realigned variants
        raVariant.setVariantDetails(VariantData.NO_LOCAL_PHASE_SET, variant.microhomology(), variant.repeatSequence(), variant.repeatCount());

        return raVariant;
    }

    public static int codonForBase(int codingBase)
    {
        // both coding base and amino acid position start at 1, so eg coding base of 1-3 = amino acid 1, 4-6 = 2 etc
        // this CodonIndex may precede the first alt-base for INDELs for the reason described above
        return  (codingBase - 1) / 3 + 1;
    }
}
