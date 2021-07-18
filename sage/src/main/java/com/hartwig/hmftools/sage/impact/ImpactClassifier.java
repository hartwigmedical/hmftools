package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.codon.Codons.START_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.START_CODON;
import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.isStopCodon;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.getCodingBaseRanges;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.FIVE_PRIME_UTR_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.THREE_PRIME_UTR_EFFECT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.FRAMESHIFT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INTRON_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.MISSENSE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_REGION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.START_LOST;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_LOST;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SYNONYMOUS_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.UPSTREAM_GENE_VARIANT;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.GENE_UPSTREAM_DISTANCE;
import static com.hartwig.hmftools.sage.impact.SpliceClassifier.isWithinSpliceRegion;

import java.util.List;

import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.CodingBaseData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.VariantConsequence;

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
        // non-coding transcripts are ignored
        if(transData.CodingStart == null)
            return null;

        int position = variant.Position;

        if(isOutsideTransRange(transData, position))
            return null;

        VariantTransImpact transImpact = checkNonCodingImpact(variant, transData);

        if(transImpact != null)
            return transImpact;

        boolean inSpliceRegion = false;

        // check intron variant
        int exonCount = transData.exons().size();
        for(int i = 0; i < exonCount; ++i)
        {
            ExonData exon = transData.exons().get(i);

            ExonData nextExon = i < exonCount - 1 ? transData.exons().get(i + 1) : null;

            if(isWithinSpliceRegion(variant, transData, exon))
            {
                inSpliceRegion = true;
                transImpact = mSpliceClassifier.classifyVariant(variant, transData, exon);

                if(transImpact != null)
                    break;
            }
            else if(nextExon != null && isWithinSpliceRegion(variant, transData, nextExon))
            {
                inSpliceRegion = true;
                transImpact = mSpliceClassifier.classifyVariant(variant, transData, nextExon);

                if(transImpact != null)
                    break;
            }

            if(exon.Start <= position && position <= exon.End)
            {
                transImpact = classifyExonicPosition(variant, transData, exon);
                break;
            }

            if(nextExon != null && position > exon.End && position < nextExon.Start)
            {
                transImpact = new VariantTransImpact(transData, INTRON_VARIANT);
                break;
            }
        }

        if(transImpact != null)
        {
            if(inSpliceRegion)
            {
                transImpact.markSpliceRegion();
                transImpact.addConsequence(SPLICE_REGION_VARIANT.description());
            }

            checkStopStartCodons(transImpact);
        }

        return transImpact;
    }

    private boolean isOutsideTransRange(final TranscriptData transData, int position)
    {
        if(transData.posStrand())
            return position < transData.TransStart - GENE_UPSTREAM_DISTANCE || position > transData.TransEnd;
        else
            return position < transData.TransStart || position > transData.TransEnd + GENE_UPSTREAM_DISTANCE;
    }

    private VariantTransImpact checkNonCodingImpact(final VariantData variant, final TranscriptData transData)
    {
        int position = variant.Position;

        if(transData.posStrand())
        {
            // check pre-gene region
            if(position >= transData.TransStart - GENE_UPSTREAM_DISTANCE && position < transData.TransStart)
                return new VariantTransImpact(transData, UPSTREAM_GENE_VARIANT);

            // check 5' and 3' UTR
            if(position >= transData.TransStart && position < transData.CodingStart)
                return new VariantTransImpact(transData, FIVE_PRIME_UTR_EFFECT);

            if(position > transData.CodingEnd && position <= transData.TransEnd)
                return new VariantTransImpact(transData, THREE_PRIME_UTR_EFFECT);

        }
        else
        {
            if(position > transData.TransEnd && position <= transData.TransEnd + GENE_UPSTREAM_DISTANCE)
                return new VariantTransImpact(transData, UPSTREAM_GENE_VARIANT);

            if(position >= transData.TransStart && position < transData.CodingStart)
                return new VariantTransImpact(transData, THREE_PRIME_UTR_EFFECT);

            if(position > transData.CodingEnd && position <= transData.TransEnd)
                return new VariantTransImpact(transData, FIVE_PRIME_UTR_EFFECT);
        }

        return null;
    }

    private VariantTransImpact classifyExonicPosition(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        if(variant.isIndel())
        {
            // for a DEL on the reverse strand this is position + baseDiff - 1
            // for an insert on the reverse strand this is position + 1

            if((variant.baseDiff() % 3) == 0)
            {
                // could use phasing info to set conservative vs disruptive type
                if(variant.isDeletion())
                    return new VariantTransImpact(transData, INFRAME_DELETION);
                else
                    return new VariantTransImpact(transData, INFRAME_INSERTION);
            }
            else
            {
                return new VariantTransImpact(transData, FRAMESHIFT_VARIANT);
            }
        }
        else
        {
            // find the start of the codon in which the first base of the variant is in, and likewise the end of the last codon
            int varLength = variant.Ref.length();
            boolean posStrand = transData.posStrand();

            int upstreamStartPos = posStrand ? variant.Position : variant.Position + varLength - 1;
            final CodingBaseData cbData = calcCodingBases(transData, upstreamStartPos);

            // this is 0 if the variant starts on the first base of a codon (phase=1), 1 if it starts on the 2nd, and 2 if on the 3rd/last
            int openCodonBases = TranscriptUtils.getUpstreamOpenCodonBases(cbData.Phase, transData.Strand, variant.baseDiff());

            int codingBaseLen = (int)(Math.ceil((openCodonBases + varLength) / 3.0)) * 3;

            int codonStartPos = posStrand ? upstreamStartPos - openCodonBases : upstreamStartPos + openCodonBases;

            List<int[]> codingRanges = getCodingBaseRanges(transData, codonStartPos, posStrand, codingBaseLen);
            String refCodingBases = mRefGenome.getBaseString(variant.Chromosome, codingRanges);

            String altCodingBases;
            String refAminoAcids;
            String altAminoAcids;

            if(posStrand)
            {
                if(openCodonBases > 0)
                {
                    altCodingBases = refCodingBases.substring(0, openCodonBases) + variant.Alt
                            + refCodingBases.substring(openCodonBases + varLength);
                }
                else
                {
                    altCodingBases = variant.Alt + refCodingBases.substring(varLength);
                }

                refAminoAcids = Codons.aminoAcidFromBases(refCodingBases);
                altAminoAcids = Codons.aminoAcidFromBases(altCodingBases);
            }
            else
            {
                if(openCodonBases > 0)
                {
                    altCodingBases = refCodingBases.substring(0, codingBaseLen - varLength - openCodonBases) + variant.Alt
                            + refCodingBases.substring(codingBaseLen - openCodonBases);
                }
                else
                {
                    altCodingBases = refCodingBases.substring(0, codingBaseLen - varLength) + variant.Alt;
                }

                refAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(refCodingBases));
                altAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(altCodingBases));
            }

            VariantConsequence consequence = refAminoAcids.equals(altAminoAcids) ? SYNONYMOUS_VARIANT : MISSENSE_VARIANT;
            VariantTransImpact transImpact = new VariantTransImpact(transData, consequence);
            transImpact.setAminoAcids("", refAminoAcids, altAminoAcids, "");
            return transImpact;
        }
    }

    private void checkStopStartCodons(final VariantTransImpact transImpact)
    {
        if(transImpact.wildtypeAA().isEmpty() || transImpact.novelAA().isEmpty())
            return;

        VariantConsequence ssConsequence = checkStopStartCodons(transImpact.wildtypeAA(), transImpact.novelAA());

        if(ssConsequence != null)
        {
            transImpact.addConsequence(ssConsequence.description());
        }
    }

    public static VariantConsequence checkStopStartCodons(final String refAminoAcids, final String altAminoAcids)
    {
        if(refAminoAcids.isEmpty() || altAminoAcids.isEmpty())
            return null;

        if(refAminoAcids.charAt(0) == START_AMINO_ACID && altAminoAcids.charAt(0) != START_AMINO_ACID)
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