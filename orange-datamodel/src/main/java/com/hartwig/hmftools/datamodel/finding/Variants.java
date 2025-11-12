package com.hartwig.hmftools.datamodel.finding;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;

import org.jetbrains.annotations.NotNull;

public final class Variants
{
    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    @NotNull
    public static List<SmallVariant> sort(@NotNull List<SmallVariant> variants)
    {
        return variants.stream().sorted((variant1, variant2) ->
        {
            double driverLikelihood1 = variant1.driverLikelihood() != null ? variant1.driverLikelihood() : -1;
            double driverLikelihood2 = variant2.driverLikelihood() != null ? variant2.driverLikelihood() : -1;

            int driverCompare = Double.compare(driverLikelihood2, driverLikelihood1);
            if(driverCompare != 0)
            {
                return driverCompare;
            }

            int geneCompare = variant1.gene().compareTo(variant2.gene());
            if(geneCompare != 0)
            {
                return geneCompare;
            }

            if(variant1.affectedCodon() == null && variant2.affectedCodon() == null)
            {
                return 0;
            }
            else if(variant1.affectedCodon() == null)
            {
                return 1;
            }
            else if(variant2.affectedCodon() == null)
            {
                return -1;
            }
            else
            {
                return Integer.compare(variant1.affectedCodon(), variant2.affectedCodon());
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static String variantField(@NotNull SmallVariant variant, boolean useSingleLetterAminoAcid)
    {
        String addon = "";

        if(!variant.isCanonical())
        {
            addon = " (alt)";
        }

        return variant.gene() + addon + " " + impactField(variant.transcriptImpact(), useSingleLetterAminoAcid);
    }

    @NotNull
    public static String hotspotField(@NotNull SmallVariant variant)
    {
        switch(variant.hotspot())
        {
            case HOTSPOT:
                return "Yes";
            case NEAR_HOTSPOT:
                return "Near";
            default:
                return "No";
        }
    }

    @NotNull
    public static String impactField(@NotNull PurpleTranscriptImpact impact, boolean useSingleLetterAminoAcid)
    {
        String hgvsProteinImpact = impact.hgvsProteinImpact();
        if(!hgvsProteinImpact.isEmpty() && !hgvsProteinImpact.equals("p.?"))
        {
            if(useSingleLetterAminoAcid)
            {
                hgvsProteinImpact = forceSingleLetterProteinAnnotation(hgvsProteinImpact);
            }
            return hgvsProteinImpact;
        }

        String hgvsCodingImpact = impact.hgvsCodingImpact();
        if(!hgvsCodingImpact.isEmpty())
        {
            return impact.codingEffect() == PurpleCodingEffect.SPLICE ? hgvsCodingImpact + " splice" : hgvsCodingImpact;
        }

        Set<PurpleVariantEffect> effects = impact.effects();
        if(effects.contains(PurpleVariantEffect.UPSTREAM_GENE))
        {
            return "upstream";
        }

        StringJoiner joiner = new StringJoiner(", ");
        for(PurpleVariantEffect effect : effects)
        {
            joiner.add(effect.effect());
        }
        return joiner.toString();
    }

    @NotNull
    public static String biallelicLikelihoodField(@NotNull SmallVariant variant)
    {
        return PERCENTAGE_FORMAT.format(variant.biallelicProbability() * 100);
    }

    @NotNull
    public static String driverLikelihoodField(@NotNull SmallVariant variant)
    {
        return variant.driverLikelihood() != null ? PERCENTAGE_FORMAT.format(variant.driverLikelihood() * 100) : "";
    }

    @NotNull
    public static String clonalLikelihoodField(@NotNull SmallVariant variant)
    {
        return PERCENTAGE_FORMAT.format(100 * variant.clonalLikelihood());
    }

    @NotNull
    public static String rnaDepthField(@NotNull SmallVariant variant, @NotNull String notAvailableString)
    {
        PurpleAllelicDepth rnaDepth = variant.rnaDepth();

        if(rnaDepth == null)
        {
            return notAvailableString;
        }

        String vafAddon = "";
        if(rnaDepth.totalReadCount() > 0)
        {
            double vaf = rnaDepth.alleleReadCount() / (double) rnaDepth.totalReadCount();
            vafAddon = " (" + PERCENTAGE_FORMAT.format(vaf * 100) + ")";
        }

        return rnaDepth.alleleReadCount() + "/" + rnaDepth.totalReadCount() + vafAddon;
    }

    @NotNull
    public static String phaseSetField(@NotNull SmallVariant variant)
    {
        List<Integer> localPhaseSets = variant.localPhaseSets();
        if(localPhaseSets == null || localPhaseSets.isEmpty())
        {
            return "";
        }

        StringJoiner joiner = new StringJoiner(", ");
        for(Integer localPhaseSet : localPhaseSets)
        {
            joiner.add(String.valueOf(localPhaseSet));
        }
        return joiner.toString();
    }

    private static String forceSingleLetterProteinAnnotation(final String proteinAnnotation)
    {
        final List<Map.Entry<String, String>> TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER = new ArrayList<>();
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Ala", "A")); // Alanine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Cys", "C")); // Cysteine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Asp", "D")); // Aspartic Acid
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Glu", "E")); // Glutamic Acid
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Phe", "F")); // Phenylalanine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Gly", "G")); // Glycine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("His", "H")); // Histidine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Ile", "I")); // Isoleucine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Lys", "K")); // Lysine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Leu", "L")); // Leucine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Met", "M")); // Methionine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Asn", "N")); // Asparagine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Pro", "P")); // Proline
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Gln", "Q")); // Glutamine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Arg", "R")); // Arginine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Ser", "S")); // Serine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Thr", "T")); // Threonine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Val", "V")); // Valine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Trp", "W")); // Tryptophan
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.add(Map.entry("Tyr", "Y")); // Tyrosine

        String convertedProteinAnnotation = proteinAnnotation;
        for(Map.Entry<String, String> mapping : TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER)
        {
            convertedProteinAnnotation = convertedProteinAnnotation.replaceAll(mapping.getKey(), mapping.getValue());
        }
        return convertedProteinAnnotation;
    }
}
