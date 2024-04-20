package com.hartwig.hmftools.pave.impact;

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
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.isSplice;
import static com.hartwig.hmftools.pave.impact.HgvsProtein.HGVS_SPLICE_UNKNOWN;
import static com.hartwig.hmftools.pave.impact.HgvsProtein.reportProteinImpact;
import static com.hartwig.hmftools.pave.impact.PaveUtils.withinTransRange;
import static com.hartwig.hmftools.pave.impact.SpliceClassifier.checkStraddlesSpliceRegion;
import static com.hartwig.hmftools.pave.impact.SpliceClassifier.isInsertIntoExonStart;
import static com.hartwig.hmftools.pave.impact.SpliceClassifier.isWithinSpliceRegion;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.pave.VariantData;

public class ImpactClassifier
{
    private final RefGenomeInterface mRefGenome;
    private final PhasedVariantClassifier mPhasedVariants;

    public ImpactClassifier(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
        mPhasedVariants = new PhasedVariantClassifier();
    }

    public PhasedVariantClassifier phasedVariants() { return mPhasedVariants; }
    public RefGenomeInterface refGenome() { return mRefGenome; }

    public VariantTransImpact classifyVariant(final VariantData variant, final TranscriptData transData)
    {
        mPhasedVariants.checkAddVariant(variant);

        if(!withinTransRange(transData, variant.Position, variant.EndPosition)) // a conservative check by 1 base for INDEL at boundary
            return null;

        VariantTransImpact transImpact = new VariantTransImpact(transData);

        CodingUtils.determineContext(variant, transData, transImpact.codingContext());
        HgvsCoding.set(variant, transImpact.codingContext());

        if(transImpact.codingContext().isCoding())
        {
            ProteinContext proteinContext = ProteinUtils.determineContext(variant, transImpact.codingContext(), transData, mRefGenome);
            transImpact.setProteinContext(proteinContext);
        }

        checkPrePostCodingImpact(variant, transImpact); // was previously a return but need to check for other effects

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
                SpliceClassifier.classifyVariant(variant, transImpact, exon, mRefGenome);
            }
            else if(!inSpliceRegion && nextExon != null && isWithinSpliceRegion(variant, transData, nextExon))
            {
                inSpliceRegion = true;
                SpliceClassifier.classifyVariant(variant, transImpact, nextExon, mRefGenome);
            }

            if(variant.altPositionsOverlap(exon.Start, exon.End) || isInsertIntoExonStart(variant, exon, transData.posStrand()))
            {
                classifyExonicPosition(variant, transImpact);
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

        if(transImpact.hasProteinContext())
        {
            List<VariantEffect> proteinEffects = transImpact.effects().stream().filter(x -> reportProteinImpact(x)).collect(Collectors.toList());

            if(!proteinEffects.isEmpty())
                transImpact.proteinContext().Hgvs = HgvsProtein.generate(transImpact.proteinContext(), proteinEffects);
        }

        // override for splice-effect variants
        if(isSplice(transImpact.topEffect()) && !transData.nonCoding())
        {
            if(!transImpact.hasProteinContext())
                transImpact.setProteinContext(new ProteinContext());

            transImpact.proteinContext().Hgvs = HGVS_SPLICE_UNKNOWN;
        }

        return transImpact;
    }

    private boolean checkPrePostCodingImpact(final VariantData variant, final VariantTransImpact transImpact)
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

    private void classifyExonicPosition(final VariantData variant, final VariantTransImpact transImpact)
    {
        final TranscriptData transData = transImpact.TransData;

        if(transData.nonCoding())
        {
            transImpact.addEffect(NON_CODING_TRANSCRIPT);
            return;
        }

        if(transImpact.hasEffect(FIVE_PRIME_UTR) || transImpact.hasEffect(THREE_PRIME_UTR))
            return;

        if(variant.isIndel())
        {
            if(transImpact.codingContext().IsFrameShift)
            {
                // could use phasing info to set conservative vs disruptive type
                transImpact.addEffect(FRAMESHIFT);
            }
            else
            {
                if(variant.isDeletion())
                    transImpact.addEffect(INFRAME_DELETION);
                else
                    transImpact.addEffect(INFRAME_INSERTION);
            }
        }
        else
        {
            if(transImpact.hasAminoAcids() && transImpact.proteinContext().hasProteinChange())
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
        if(!transImpact.hasAminoAcids())
            return;

        VariantEffect ssEffect = checkStopStartCodons(
                transImpact.proteinContext().CodonIndex, transImpact.proteinContext().RefAminoAcids, transImpact.proteinContext().AltAminoAcids);

        if(ssEffect != null)
        {
            // remove any superceded effects
            if(ssEffect == STOP_LOST || ssEffect == START_LOST)
            {
                transImpact.effects().remove(MISSENSE);
                transImpact.effects().remove(FRAMESHIFT);
            }

            if(ssEffect == STOP_GAINED)
                transImpact.effects().remove(MISSENSE);

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

    public static VariantTransImpact selectAlignedImpacts(final VariantTransImpact transImpact, final VariantTransImpact raTransImpact)
    {
        if(raTransImpact == null)
            return transImpact;

        // special exceptions - cannot attempt right realignment if a stop-codon was created
        if(transImpact.hasEffect(STOP_GAINED) && transImpact.TransData.posStrand())
            return transImpact;

        // take the least impact effect but favour non-splice over splice
        boolean nonRaHasSplice = transImpact.effects().stream().anyMatch(x -> isSplice(x));
        boolean raHasSplice = raTransImpact.effects().stream().anyMatch(x -> isSplice(x));

        boolean useRealigned = false;

        if(raHasSplice || nonRaHasSplice)
        {
            if(nonRaHasSplice && !raHasSplice)
                useRealigned = true;
        }
        else if(raTransImpact.topRank() != transImpact.topRank())
        {
            useRealigned = raTransImpact.topRank() < transImpact.topRank();
        }
        else if(raTransImpact.proteinContext() != null && transImpact.proteinContext() != null
        && raTransImpact.proteinContext().IsDuplication != transImpact.proteinContext().IsDuplication)
        {
            useRealigned = raTransImpact.proteinContext().IsDuplication;
        }

        if(useRealigned)
        {
            raTransImpact.markRealigned();

            if(transImpact.spliceImpactType() != SpliceImpactType.UNKNOWN && raTransImpact.spliceImpactType() == SpliceImpactType.UNKNOWN)
                raTransImpact.setSpliceImpactType(transImpact.spliceImpactType());

            return raTransImpact;
        }
        else
        {
            return transImpact;
        }
    }

    public List<VariantData> processPhasedVariants(int currentLocalPhaseSet)
    {
        if(!mPhasedVariants.hasCompleteVariants(currentLocalPhaseSet))
            return null;

        List<PhasedVariants> phasedVariantList = mPhasedVariants.popCompletePhasedVariants(currentLocalPhaseSet);
        List<VariantData> variants = Lists.newArrayList();

        for(PhasedVariants phasedVariants : phasedVariantList)
        {
            mPhasedVariants.reclassifyPhasedVariants(phasedVariants, mRefGenome);
            variants.addAll(phasedVariants.variants());
        }

        return variants;
    }

}