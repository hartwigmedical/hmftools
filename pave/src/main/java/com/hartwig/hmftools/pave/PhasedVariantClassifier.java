package com.hartwig.hmftools.pave;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.START_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.pave.HgvsProtein.reportProteinImpact;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.ProteinUtils.trimAminoAcids;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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

        // PV_LOGGER.trace("processing {} phased variant variants with localPhaseId({})", variants.Variants.size(), variants.LocalPhaseId);

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

        if(transImpacts.stream().anyMatch(x -> x.hasEffect(STOP_LOST) || x.hasEffect(START_LOST)))
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
            PV_LOGGER.trace("lps({}) varCount({}) out-of-frame with indel base count({})",
                    localPhaseSet, variants.size(), indelBaseTotal);
            return;
        }

        final String chromosome = variants.get(0).Chromosome;

        String combinedRefCodons = "";
        String combinedAltCodons = "";
        int minCodonIndex = -1;
        int lastRefCodonEnd = 0;
        boolean hasOverlaps = false;

        for(int i = 0; i < transImpacts.size(); ++i)
        {
            VariantData variant = variants.get(i);
            VariantTransImpact transImpact = transImpacts.get(i);

            if(!transImpact.hasProteinContext())
                continue;

            // ignore SNVs/MNVs if they don't overlap the same codons as any of the INDELs
            int refCodonStart = transImpact.proteinContext().refCodingBaseStart();
            int refCodonEnd = transImpact.proteinContext().refCodingBaseEnd();

            VariantTransImpact nextTransImpact = i < transImpacts.size() - 1 ? transImpacts.get(i + 1) : null;

            int nextRefCodonStart = nextTransImpact != null && nextTransImpact.hasProteinContext()
                    ? nextTransImpact.proteinContext().refCodingBaseStart() : 0;

            boolean overlapsOnStart = lastRefCodonEnd > 0 && refCodonStart <= lastRefCodonEnd;
            boolean overlapsOnEnd = nextRefCodonStart > 0 && refCodonEnd >= nextRefCodonStart;

            hasOverlaps = overlapsOnStart || overlapsOnEnd;

            if(variant.isBaseChange() && !overlapsOnStart && !overlapsOnEnd)
            {
                PV_LOGGER.trace("lps({}) var({}) skip no overlap range({} -> {}) vs last({}) and next({})",
                        localPhaseSet, variant, refCodonStart, refCodonEnd, lastRefCodonEnd, nextRefCodonStart);
                continue;
            }

            PV_LOGGER.trace("lps({}) var({}) codons({} -> {}) range({} - {}) overlaps(start={} end={})",
                    localPhaseSet, variant, transImpact.proteinContext().RefCodonBases, transImpact.proteinContext().AltCodonBases,
                    refCodonStart, refCodonEnd, overlapsOnStart, overlapsOnEnd);

            VariantData prevVariant = i > 0 ? variants.get(i - 1) : null;

            minCodonIndex = i == 0 ?
                    transImpact.proteinContext().CodonIndex : min(transImpact.proteinContext().CodonIndex, minCodonIndex);

            if(overlapsOnStart && lastRefCodonEnd == refCodonEnd && variant.isBaseChange() && !prevVariant.isBaseChange()
            && transImpact.proteinContext().RefCodonBases.length() == 3) // at most one codon for now
            {
                // no need to re-add the ref bases
                // take an modified bases that are after the end of the previous variant
                int posPhase = transImpact.codingContext().UpstreamPhase;
                if(posPhase == PHASE_0)
                    posPhase = 3;

                // phase now goes 1, 2, 3
                String lastAltCodon = combinedAltCodons.substring(combinedAltCodons.length() - 3);
                String altCodon = transImpact.proteinContext().AltCodonBases;
                String mixedAltCodon = "";

                for(int j = 0; j < 3; ++j)
                {
                    if(j + 1 < posPhase)
                        mixedAltCodon += lastAltCodon.charAt(j);
                    else
                        mixedAltCodon += altCodon.charAt(j);
                }

                combinedAltCodons = combinedAltCodons.substring(0, combinedAltCodons.length() - 3) + mixedAltCodon;
                hasOverlaps = false;
            }
            else
            {
                combinedRefCodons += transImpact.proteinContext().RefCodonBases;

                // fill in any missing gaps in the codons
                if(lastRefCodonEnd > 0 && !overlapsOnStart && refCodonStart > lastRefCodonEnd + 1)
                {
                    String gapCodonBases = refGenome.getBaseString(chromosome, lastRefCodonEnd + 1, refCodonStart - 1);
                    combinedRefCodons += gapCodonBases;
                    combinedAltCodons += gapCodonBases;
                }

                combinedAltCodons += transImpact.proteinContext().AltCodonBases;
            }

            lastRefCodonEnd = transImpact.proteinContext().refCodingBaseEnd();
        }

        if(hasOverlaps && indelBaseTotal == 0)
            return;

        if(!isCodonMultiple(combinedRefCodons.length()) || !isCodonMultiple(combinedAltCodons.length()))
            return;

        ProteinContext combinedPc = new ProteinContext();
        combinedPc.CodonIndex = minCodonIndex;
        combinedPc.RefCodonBases = combinedRefCodons;
        combinedPc.AltCodonBases = combinedAltCodons;

        // get any additional downstream bases to complete the alt codon(s)?

        // check synonymous vs missense from combined bases
        final TranscriptData transData = transImpacts.get(0).TransData;

        VariantEffect combinedEffect;

        if(transData.posStrand())
        {
            combinedPc.RefAminoAcids = Codons.aminoAcidFromBases(combinedRefCodons);
            combinedPc.AltAminoAcids = Codons.aminoAcidFromBases(combinedAltCodons);
        }
        else
        {
            combinedPc.RefAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(combinedRefCodons));
            combinedPc.AltAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(combinedAltCodons));
        }

        trimAminoAcids(combinedPc, true, true, true);

        if(indelBaseTotal == 0)
        {
            combinedEffect = combinedPc.RefAminoAcids.equals(combinedPc.AltAminoAcids) ? SYNONYMOUS : MISSENSE;

            PV_LOGGER.trace("lps({}) varCount({}) combinedEffect({}) from aminoAcids({} -> {})",
                    localPhaseSet, variants.size(), combinedEffect, combinedPc.RefAminoAcids, combinedPc.AltAminoAcids);
        }
        else
        {
            combinedEffect = indelBaseTotal > 0 ? PHASED_INFRAME_INSERTION : PHASED_INFRAME_DELETION;

            PV_LOGGER.trace("lps({}) varCount({}) indelBaseTotal({})",
                    localPhaseSet, variants.size(), indelBaseTotal);
        }

        combinedPc.Hgvs = HgvsProtein.generate(combinedPc, combinedEffect);


        for(VariantTransImpact transImpact : transImpacts)
        {
            if(!transImpact.hasProteinContext())
                continue;

            transImpact.markPhasedFrameshift();
            transImpact.setProteinContext(combinedPc);

            transImpact.codingContext().IsFrameShift = false;
            transImpact.effects().remove(FRAMESHIFT);
            transImpact.addEffect(combinedEffect);
        }
    }
}
