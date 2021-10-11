package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

public class PhasedVariantClassifier
{
    private final List<PhasedVariants> mPhasedVariants;

    public PhasedVariantClassifier()
    {
        mPhasedVariants = Lists.newArrayList();
    }

    public void checkAddVariant(final VariantData variant)
    {
        if(!variant.hasLocalPhaseSet())
            return;

        PhasedVariants phasedVariants = mPhasedVariants.stream()
                .filter(x -> x.LocalPhaseId == variant.localPhaseSet()).findFirst().orElse(null);

        if(phasedVariants == null)
        {
            phasedVariants = new PhasedVariants(variant.localPhaseSet());
            mPhasedVariants.add(phasedVariants);
        }

        phasedVariants.addVariant(variant);
    }

    public boolean hasCompleteVariants(int currentLocalPhaseSet)
    {
        return !mPhasedVariants.isEmpty() && mPhasedVariants.get(0).LocalPhaseId != currentLocalPhaseSet;
    }

    public List<PhasedVariants> popCompletePhasedVariants(int currentLocalPhaseSet)
    {
        List<PhasedVariants> completeVariants = Lists.newArrayList();

        while(!mPhasedVariants.isEmpty())
        {
            if(mPhasedVariants.get(0).LocalPhaseId != currentLocalPhaseSet)
            {
                completeVariants.add(mPhasedVariants.get(0));
                mPhasedVariants.remove(0);
            }
            else
            {
                break;
            }
        }

        return completeVariants;
    }

    public List<PhasedVariants> allPhasedVariants() { return mPhasedVariants; }
    public void clear() { mPhasedVariants.clear(); }

    public static void reclassifyPhasedVariants(final PhasedVariants phasedVariants, final RefGenomeInterface refGenome)
    {
        if(phasedVariants.variants().size() < 2)
            return;

        // ignore groups without INDELS
        if(phasedVariants.variants().stream().noneMatch(x -> x.isIndel()))
            return;

        // PV_LOGGER.debug("processing {} phased variant variants with localPhaseId({})", variants.Variants.size(), variants.LocalPhaseId);

        // assume that variants are ordered  by ascending position
        VariantData variant = phasedVariants.variants().get(0);

        if(variant.getImpacts().isEmpty())
            return;

        // process impacts by matching gene and transcript
        for(Map.Entry<String,List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
        {
            String gene = entry.getKey();

            for(VariantTransImpact transImpact : entry.getValue())
            {
                List<VariantTransImpact> transImpacts = Lists.newArrayList(transImpact);
                List<VariantData> variants = Lists.newArrayList(variant);

                for(int i = 1; i < phasedVariants.variants().size(); ++i)
                {
                    VariantData nextVariant = phasedVariants.variants().get(i);

                    List<VariantTransImpact> otherImpacts = nextVariant.getImpacts().get(gene);

                    if(otherImpacts != null)
                    {
                        VariantTransImpact nextImpact = otherImpacts.stream()
                                .filter(x -> x.TransData.TransName.equals(transImpact.TransData.TransName))
                                .findFirst().orElse(null);

                        if(nextImpact != null)
                        {
                            transImpacts.add(nextImpact);
                            variants.add(nextVariant);
                        }
                    }
                }

                reclassifyImpacts(phasedVariants.LocalPhaseId, variants, transImpacts, refGenome);
            }
        }
    }

    private static void reclassifyImpacts(
            int localPhaseSet, final List<VariantData> variants, final List<VariantTransImpact> transImpacts, final RefGenomeInterface refGenome)
    {
        // ignore if not 2 or more with coding impacts
        if(transImpacts.stream().filter(x -> x.hasProteinContext()).count() < 2)
            return;

        // ignore if not phased anyway
        int indelBaseTotal = 0;

        for(int i = 0; i < transImpacts.size(); ++i)
        {
            VariantTransImpact transImpact = transImpacts.get(i);
            VariantData variant = variants.get(i);

            if(!transImpact.hasProteinContext() || !variant.isIndel())
                continue;

            indelBaseTotal += variant.isInsert() ? variant.baseDiff() : -transImpact.codingContext().DeletedCodingBases;
        }

        if(!isCodonMultiple(indelBaseTotal))
        {
            PV_LOGGER.debug("lps({}) varCount({}) out-of-frame with indel base count({})",
                    localPhaseSet, variants.size(), indelBaseTotal);
            return;
        }

        final String chromosome = variants.get(0).Chromosome;

        String combinedRefCodons = "";
        String combinedAltCodons = "";
        int lastRefCodonEnd = 0;
        boolean hasOverlaps = false;

        for(int i = 0; i < transImpacts.size(); ++i)
        {
            VariantData variant = variants.get(i);
            VariantTransImpact transImpact = transImpacts.get(i);

            if(!transImpact.hasProteinContext())
                continue;

            // ignore SNVs/MNVs if they don't overlap the same codons as any of the INDELs
            int refCodonStart = transImpact.proteinContext().RefCodonsRange[SE_START];
            int refCodonEnd = transImpact.proteinContext().RefCodonsRange[SE_END];

            VariantTransImpact nextTransImpact = i < transImpacts.size() - 1 ? transImpacts.get(i + 1) : null;

            int nextRefCodonStart = nextTransImpact != null && nextTransImpact.hasProteinContext()
                    ? nextTransImpact.proteinContext().RefCodonsRange[SE_START] : 0;

            boolean overlapsOnStart = lastRefCodonEnd > 0 && refCodonStart <= lastRefCodonEnd;
            boolean overlapsOnEnd = nextRefCodonStart > 0 && refCodonStart >= nextRefCodonStart;

            hasOverlaps = overlapsOnStart || overlapsOnEnd;

            if(variant.isBaseChange() && !overlapsOnStart && !overlapsOnEnd)
            {
                PV_LOGGER.debug("lps({}) var({}) skip no overlap range({} -> {}) vs last({}) and next({})",
                        localPhaseSet, variant, refCodonStart, refCodonEnd, lastRefCodonEnd, nextRefCodonStart);
                continue;
            }

            combinedRefCodons += transImpact.proteinContext().RefCodonBases;

            PV_LOGGER.debug("lps({}) var({}) codons({} -> {}) range({} - {}) overlaps(start={} end={})",
                    localPhaseSet, variant, transImpact.proteinContext().RefCodonBases, transImpact.proteinContext().AltCodonBases,
                    refCodonStart, refCodonEnd, overlapsOnStart, overlapsOnEnd);

            // fill in any missing gaps in the codons
            if(lastRefCodonEnd > 0 && !overlapsOnStart && refCodonStart > lastRefCodonEnd + 1)
            {
                String gapCodonBases = refGenome.getBaseString(chromosome, lastRefCodonEnd + 1, refCodonStart - 1);
                combinedRefCodons += gapCodonBases;
                combinedAltCodons += gapCodonBases;
            }

            lastRefCodonEnd = transImpact.proteinContext().RefCodonsRange[SE_END];
            combinedAltCodons += transImpact.proteinContext().AltCodonBases;
        }

        if(hasOverlaps)
            return;

        // check synonymous vs missense from combined bases
        final TranscriptData transData = transImpacts.get(0).TransData;

        String refAminoAcids = "";
        String altAminoAcids = "";

        if(transData.posStrand())
        {
            refAminoAcids = Codons.aminoAcidFromBases(combinedRefCodons);
            altAminoAcids = Codons.aminoAcidFromBases(combinedAltCodons);
        }
        else
        {
            refAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(combinedRefCodons));
            altAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(combinedAltCodons));
        }

        VariantEffect combinedEffect;

        if(indelBaseTotal == 0)
        {
            combinedEffect = refAminoAcids.equals(altAminoAcids) ? SYNONYMOUS : MISSENSE;
        }
        else
        {
            combinedEffect = indelBaseTotal > 0 ? INFRAME_INSERTION : INFRAME_DELETION;
        }

        PV_LOGGER.debug("lps({}) varCount({}) combinedEffect({}) from aminoAcids({} -> {})",
                localPhaseSet, variants.size(), combinedEffect, refAminoAcids, altAminoAcids);

        for(VariantTransImpact transImpact : transImpacts)
        {
            if(!transImpact.hasProteinContext())
                continue;

            transImpact.proteinContext().RefAminoAcids = refAminoAcids;
            transImpact.proteinContext().AltAminoAcids = altAminoAcids;
            transImpact.codingContext().IsFrameShift = false;
            transImpact.effects().remove(FRAMESHIFT);
            transImpact.addEffect(combinedEffect);
        }
    }
}
