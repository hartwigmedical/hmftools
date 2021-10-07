package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.START_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FIVE_PRIME_UTR;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INTRONIC;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.NON_CODING_TRANSCRIPT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.START_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.THREE_PRIME_UTR;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.UPSTREAM_GENE;
import static com.hartwig.hmftools.pave.CodingUtils.determineContext;
import static com.hartwig.hmftools.pave.PaveUtils.withinTransRange;
import static com.hartwig.hmftools.pave.SpliceClassifier.checkStraddlesSpliceRegion;
import static com.hartwig.hmftools.pave.SpliceClassifier.isWithinSpliceRegion;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

public class ImpactClassifier
{
    private final RefGenomeInterface mRefGenome;
    private final SpliceClassifier mSpliceClassifier;

    public ImpactClassifier(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
        mSpliceClassifier = new SpliceClassifier();
    }

    public VariantTransImpact classifyVariant(final VariantData variant, final TranscriptData transData)
    {
        if(!withinTransRange(transData, variant.Position, variant.EndPosition)) // a conservative check by 1 base for INDEL at boundary
            return null;

        VariantTransImpact transImpact = new VariantTransImpact(transData);

        CodingContext codingContext = CodingUtils.determineContext(variant, transData);
        transImpact.setCodingContext(codingContext);

        if(codingContext.isCoding())
        {
            ProteinContext proteinContext = ProteinUtils.determineContext(variant, codingContext, transData, mRefGenome);
            transImpact.setProteinContext(proteinContext);
        }

        if(checPrePostCodingImpact(variant, transImpact))
            return transImpact;

        boolean inSpliceRegion = false;

        // loop through all exons checking if the variant falls in a splice region and is intronic or exonic
        int exonCount = transData.exons().size();
        for(int i = 0; i < exonCount; ++i)
        {
            ExonData exon = transData.exons().get(i);
            ExonData nextExon = i < exonCount - 1 ? transData.exons().get(i + 1) : null;

            if(!inSpliceRegion && isWithinSpliceRegion(variant, transData, exon))
            {
                inSpliceRegion = true;
                transImpact.addEffect(mSpliceClassifier.classifyVariant(variant, transData, exon));
            }
            else if(!inSpliceRegion && nextExon != null && isWithinSpliceRegion(variant, transData, nextExon))
            {
                inSpliceRegion = true;
                transImpact.addEffect(mSpliceClassifier.classifyVariant(variant, transData, nextExon));
            }

            if(variant.altPositionsOverlap(exon.Start, exon.End))
            {
                classifyExonicPosition(variant, transImpact, exon);
                break;
            }

            if(nextExon != null && variant.altPositionsWithin(exon.End + 1, nextExon.Start - 1))
            {
                transImpact.addEffect(INTRONIC);
                break;
            }

            // variant which full straddles the splice region
            if(!inSpliceRegion)
            {
                VariantEffect effect = checkStraddlesSpliceRegion(variant, transData, exon);
                if(effect != null)
                {
                    transImpact.addEffect(effect);
                    break;
                }
            }
        }

        if(inSpliceRegion)
            transImpact.markSpliceRegion();

        checkStopStartCodons(transImpact);

        return transImpact;
    }

    private boolean checPrePostCodingImpact(final VariantData variant, final VariantTransImpact transImpact)
    {
        final TranscriptData transData = transImpact.TransData;

        // check pre-gene region
        if((transData.posStrand() && variant.altBasesBelow(transData.TransStart))
        || (!transData.posStrand() && variant.altBasesAbove(transData.TransEnd)))
        {
            transImpact.addEffect(UPSTREAM_GENE);
            return true;
        }

        if(!transImpact.isExonic() || transData.nonCoding())
            return false;

        if(transData.posStrand())
        {
            // check 5' and 3' UTR
            if(variant.altBasesBelow(transData.CodingStart))
            {
                transImpact.addEffect(FIVE_PRIME_UTR);
                return true;
            }

            if(variant.altBasesAbove(transData.CodingEnd))
            {
                transImpact.addEffect(THREE_PRIME_UTR);
                return true;
            }
        }
        else
        {
            if(variant.altBasesBelow(transData.CodingStart))
            {
                transImpact.addEffect(THREE_PRIME_UTR);
                return true;
            }

            if(variant.altBasesAbove(transData.CodingEnd))
            {
                transImpact.addEffect(FIVE_PRIME_UTR);
                return true;
            }
        }

        return false;
    }

    private void classifyExonicPosition(final VariantData variant, final VariantTransImpact transImpact, final ExonData exon)
    {
        final TranscriptData transData = transImpact.TransData;

        if(transData.nonCoding())
        {
            transImpact.addEffect(NON_CODING_TRANSCRIPT);
            return;
        }

        if(variant.isIndel())
        {
            // if only the end of an exon is affected (either strand) then this doesn't change the coding bases
            if(variant.EndPosition == exon.Start || variant.Position == exon.End)
                return;

            // if the variant crosses the exon boundary then only check the bases which are exonic to decide between frameshift or inframe
            int exonicBases = variant.isDeletion() ? transImpact.codingContext().DeletedCodingBases : variant.baseDiff();

            if((exonicBases % 3) == 0)
            {
                // could use phasing info to set conservative vs disruptive type

                if(variant.isDeletion())
                    transImpact.addEffect(INFRAME_DELETION);
                else
                    transImpact.addEffect(INFRAME_INSERTION);
            }
            else
            {
                transImpact.addEffect(FRAMESHIFT);
            }
        }
        else
        {
            if(transImpact.hasCodingData() && transImpact.proteinContext().hasProteinChange())
            {
                if(!transImpact.proteinContext().AltAminoAcids.equals(STOP_AMINO_ACID)) // stop-gained is a special case
                    transImpact.addEffect(MISSENSE);
            }
            else
            {
                transImpact.addEffect(SYNONYMOUS);
            }
        }
    }

    private void checkStopStartCodons(final VariantTransImpact transImpact)
    {
        if(!transImpact.hasCodingData())
            return;

        VariantEffect ssEffect = checkStopStartCodons(
                transImpact.proteinContext().StartPosition, transImpact.proteinContext().RefAminoAcids, transImpact.proteinContext().AltAminoAcids);

        if(ssEffect != null)
        {
            transImpact.effects().remove(MISSENSE); // superceded if present so remove
            transImpact.addEffect(ssEffect);
        }
    }

    public static VariantEffect checkStopStartCodons(int aminoAcidStartPos, final String refAminoAcids, final String altAminoAcids)
    {
        if(refAminoAcids.isEmpty() || altAminoAcids.isEmpty())
            return null;

        if(aminoAcidStartPos == 1 && refAminoAcids.charAt(0) == START_AMINO_ACID && altAminoAcids.charAt(0) != START_AMINO_ACID)
            return START_LOST;

        if(refAminoAcids.charAt(refAminoAcids.length() - 1) != STOP_AMINO_ACID)
        {
            // has the alt gained a stop codon anywhere
            for(int i = 0; i < altAminoAcids.length(); ++i)
            {
                if(altAminoAcids.charAt(i) == STOP_AMINO_ACID)
                    return STOP_GAINED;
            }
        }
        else
        {
            boolean hasStopCodon = false;
            for(int i = 0; i < altAminoAcids.length(); ++i)
            {
                if(altAminoAcids.charAt(i) == STOP_AMINO_ACID)
                {
                    hasStopCodon = true;
                    break;
                }
            }

            if(!hasStopCodon)
                return STOP_LOST;
        }

        return null;
    }
}