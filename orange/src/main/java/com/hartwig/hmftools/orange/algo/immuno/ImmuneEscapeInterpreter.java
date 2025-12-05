package com.hartwig.hmftools.orange.algo.immuno;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.immuno.ImmutableImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ImmuneEscapeInterpreter
{
    private static final Logger LOGGER = LogManager.getLogger(ImmuneEscapeInterpreter.class);

    static final Set<String> HLA_GENES = Sets.newHashSet("HLA-A", "HLA-B", "HLA-C");
    static final Set<String> ANTIGEN_PRESENTATION_GENES = Sets.newHashSet("B2M", "CALR", "TAP1", "TAP2", "TABBP", "NLRC5", "CIITA", "RFX5");
    static final Set<String> IFN_GAMMA_PATHWAY_GENES = Sets.newHashSet("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1");
    static final Set<String> PD_L1_GENES = Sets.newHashSet("CD274");
    static final Set<String> CD58_GENES = Sets.newHashSet("CD58");
    static final Set<String> EPIGENETIC_SETDB1_GENES = Sets.newHashSet("SETDB1");

    @NotNull
    public static ImmuneEscapeRecord interpret(final PurpleRecord purple, final LinxRecord linx)
    {
        return ImmutableImmuneEscapeRecord.builder()
                .hasHlaEscape(anyGeneWithLOH(purple, HLA_GENES) || anyGeneWithInactivation(purple, linx, HLA_GENES))
                .hasAntigenPresentationPathwayEscape(anyGeneWithInactivation(purple, linx, ANTIGEN_PRESENTATION_GENES))
                .hasIFNGammaPathwayEscape(anyGeneWithInactivation(purple, linx, IFN_GAMMA_PATHWAY_GENES))
                .hasPDL1OverexpressionEscape(anyGeneWithAmplification(purple, PD_L1_GENES))
                .hasCD58InactivationEscape(anyGeneWithInactivation(purple, linx, CD58_GENES))
                .hasEpigeneticSETDB1Escape(anyGeneWithAmplification(purple, EPIGENETIC_SETDB1_GENES))
                .build();
    }

    private static boolean anyGeneWithLOH(final PurpleRecord purple, final Set<String> genesToCheck)
    {
        for(String geneToCheck : genesToCheck)
        {
            if(hasLOH(purple.somaticGeneCopyNumbers(), geneToCheck))
            {
                return true;
            }
        }

        return false;
    }

    private static boolean hasLOH(final List<PurpleGeneCopyNumber> allSomaticGeneCopyNumbers, final String geneToCheck)
    {
        for(PurpleGeneCopyNumber somaticGeneCopyNumber : allSomaticGeneCopyNumbers)
        {
            if(somaticGeneCopyNumber.gene().equals(geneToCheck))
            {
                return somaticGeneCopyNumber.minCopyNumber() > 0.5  && somaticGeneCopyNumber.minMinorAlleleCopyNumber() < 0.3;
            }
        }

        LOGGER.warn("Could not find gene copy number data for gene: {}", geneToCheck);
        return false;
    }

    private static boolean anyGeneWithInactivation(final PurpleRecord purple, final LinxRecord linx,
            final Set<String> genesToCheck)
    {
        for(String geneToCheck : genesToCheck)
        {
            boolean hasInactivationVariant = hasAnyInactivationVariant(purple.otherSomaticVariants(), geneToCheck);
            boolean hasGeneDeletion = isDeleted(purple.driverSomaticGainsDels(), geneToCheck);
            boolean hasHomozygousDisruption = isHomozygouslyDisrupted(linx.somaticHomozygousDisruptions(), geneToCheck);

            if(hasInactivationVariant || hasGeneDeletion || hasHomozygousDisruption)
            {
                return true;
            }
        }
        return false;
    }

    private static boolean hasAnyInactivationVariant(final List<PurpleVariant> allSomaticVariants, final String geneToCheck)
    {
        for(PurpleVariant somaticVariant : allSomaticVariants)
        {
            if(somaticVariant.gene().equals(geneToCheck))
            {
                PurpleCodingEffect canonicalCodingEffect = somaticVariant.canonicalImpact().codingEffect();
                boolean hasLOFImpact = canonicalCodingEffect == PurpleCodingEffect.SPLICE ||
                        canonicalCodingEffect == PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT;
                boolean hasBiallelicMissenseImpact = canonicalCodingEffect == PurpleCodingEffect.MISSENSE && somaticVariant.biallelic();

                boolean isClonal = somaticVariant.subclonalLikelihood() < 0.5;
                if(isClonal && (hasLOFImpact || hasBiallelicMissenseImpact))
                {
                    return true;
                }
            }
        }
        return false;
    }

    private static boolean isDeleted(final List<PurpleGainDeletion> allSomaticGainDels, final String geneToCheck)
    {
        for(PurpleGainDeletion somaticGainDel : allSomaticGainDels)
        {
            if(somaticGainDel.gene().equals(geneToCheck) && somaticGainDel.isCanonical())
            {
                return somaticGainDel.interpretation() == CopyNumberInterpretation.FULL_DEL
                        || somaticGainDel.interpretation() == CopyNumberInterpretation.PARTIAL_DEL;
            }
        }

        return false;
    }

    private static boolean isHomozygouslyDisrupted(
            final List<LinxHomozygousDisruption> somaticHomozygousDisruptions, final String geneToCheck)
    {
        for(LinxHomozygousDisruption somaticHomozygousDisruption : somaticHomozygousDisruptions)
        {
            if(somaticHomozygousDisruption.gene().equals(geneToCheck) && somaticHomozygousDisruption.isCanonical())
            {
                return true;
            }
        }
        return false;
    }

    private static boolean anyGeneWithAmplification(final PurpleRecord purple, final Set<String> genesToCheck)
    {
        for(String geneToCheck : genesToCheck)
        {
            if(isAmplified(purple.driverSomaticGainsDels(), geneToCheck))
            {
                return true;
            }
        }
        return false;
    }

    private static boolean isAmplified(final List<PurpleGainDeletion> somaticGainDels, final String geneToCheck)
    {
        for(PurpleGainDeletion somaticGainDel : somaticGainDels)
        {
            if(somaticGainDel.gene().equals(geneToCheck) && somaticGainDel.isCanonical())
            {
                return somaticGainDel.interpretation() == CopyNumberInterpretation.FULL_GAIN;
            }
        }
        return false;
    }
}
