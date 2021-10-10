package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

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

                reclassifyImpacts(variants, transImpacts, refGenome);
            }
        }
    }

    private static void reclassifyImpacts(
            final List<VariantData> variants, final List<VariantTransImpact> transImpacts, final RefGenomeInterface refGenome)
    {
        if(transImpacts.stream().filter(x -> x.hasProteinContext()).count() < 2)
            return;

        final String chromosome = variants.get(0).Chromosome;

        String combinedRefCodons = "";
        String combinedAltCodons = "";
        int lastCodonBase = 0;

        for(int i = 0; i < transImpacts.size(); ++i)
        {
            VariantData variant = variants.get(i);
            VariantTransImpact transImpact = transImpacts.get(i);
            combinedRefCodons += transImpact.proteinContext().RefCodonBases;

            PV_LOGGER.debug("var({}) codons({} -> {}) range({} - {}) {}",
                    variant, transImpact.proteinContext().RefCodonBases, transImpact.proteinContext().AltCodonBases,
                    transImpact.proteinContext().RefCodonsRange[SE_START], transImpact.proteinContext().RefCodonsRange[SE_END],
                    i > 0 && transImpact.proteinContext().RefCodonsRange[SE_START] <= lastCodonBase ? "overlaps" : "");

            // fill in any missing gaps in the codons
            if(i > 0)
            {
                int nextCodonBase = transImpact.proteinContext().RefCodonsRange[SE_START];

                if(nextCodonBase > lastCodonBase + 1)
                {
                    String gapCodonBases = refGenome.getBaseString(chromosome, lastCodonBase, nextCodonBase);
                    combinedRefCodons += gapCodonBases;
                    combinedAltCodons += gapCodonBases;
                }
            }

            lastCodonBase = transImpact.proteinContext().RefCodonsRange[SE_END];
            combinedAltCodons += transImpact.proteinContext().AltCodonBases;
        }
    }

}
