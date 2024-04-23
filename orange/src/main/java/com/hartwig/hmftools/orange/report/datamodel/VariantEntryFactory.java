package com.hartwig.hmftools.orange.report.datamodel;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.orange.report.interpretation.Drivers;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class VariantEntryFactory
{
    @NotNull
    public static List<VariantEntry> create(@NotNull List<PurpleVariant> variants, @NotNull List<PurpleDriver> drivers)
    {
        List<VariantEntry> entries = Lists.newArrayList();
        for(PurpleVariant variant : variants)
        {
            PurpleDriver driver = Drivers.canonicalMutationEntryForGene(drivers, variant.gene());
            entries.add(toVariantEntry(variant, driver));
        }

        for(PurpleDriver nonCanonicalDriver : Drivers.nonCanonicalMutationEntries(drivers))
        {
            List<PurpleVariant> nonCanonicalVariants = findReportedVariantsForDriver(variants, nonCanonicalDriver);
            for(PurpleVariant nonCanonicalVariant : nonCanonicalVariants)
            {
                entries.add(toVariantEntry(nonCanonicalVariant, nonCanonicalDriver));
            }
        }

        return entries;
    }

    @NotNull
    private static VariantEntry toVariantEntry(@NotNull PurpleVariant variant, @Nullable PurpleDriver driver)
    {
        PurpleTranscriptImpact transcriptImpact;

        if(driver != null)
        {
            transcriptImpact = findTranscriptImpact(variant, driver.transcript());
            if(transcriptImpact == null)
            {
                throw new IllegalStateException("Could not find impact on transcript " + driver.transcript() + " for variant " + variant);
            }
        }
        else
        {
            transcriptImpact = variant.canonicalImpact();
        }

        return ImmutableVariantEntry.builder()
                .gene(variant.gene())
                .isCanonical(driver == null || driver.transcript().equals(variant.canonicalImpact().transcript()))
                .affectedCodon(transcriptImpact.affectedCodon())
                .impact(determineImpact(transcriptImpact))
                .variantCopyNumber(variant.adjustedCopyNumber() * Math.max(0, Math.min(1, variant.adjustedVAF())))
                .totalCopyNumber(variant.adjustedCopyNumber())
                .minorAlleleCopyNumber(variant.minorAlleleCopyNumber())
                .biallelic(variant.biallelic())
                .hotspot(variant.hotspot())
                .driverLikelihood(driver != null ? driver.driverLikelihood() : null)
                .clonalLikelihood(1 - variant.subclonalLikelihood())
                .localPhaseSets(variant.localPhaseSets())
                .rnaDepth(variant.rnaDepth())
                .genotypeStatus(variant.genotypeStatus())
                .build();
    }

    @NotNull
    private static List<PurpleVariant> findReportedVariantsForDriver(@NotNull List<PurpleVariant> variants, @NotNull PurpleDriver driver)
    {
        List<PurpleVariant> reportedVariantsForDriver = Lists.newArrayList();
        List<PurpleVariant> reportedVariantsForGene = findReportedVariantsForGene(variants, driver.gene());
        for(PurpleVariant variant : reportedVariantsForGene)
        {
            if(findTranscriptImpact(variant, driver.transcript()) != null)
            {
                reportedVariantsForDriver.add(variant);
            }
        }

        return reportedVariantsForDriver;
    }

    @NotNull
    private static List<PurpleVariant> findReportedVariantsForGene(@NotNull List<PurpleVariant> variants, @NotNull String geneToFind)
    {
        List<PurpleVariant> reportedVariantsForGene = Lists.newArrayList();
        for(PurpleVariant variant : variants)
        {
            if(variant.reported() && variant.gene().equals(geneToFind))
            {
                reportedVariantsForGene.add(variant);
            }
        }
        return reportedVariantsForGene;
    }

    @Nullable
    @VisibleForTesting
    static PurpleTranscriptImpact findTranscriptImpact(@NotNull PurpleVariant variant, @NotNull String transcriptToFind)
    {
        if(variant.canonicalImpact().transcript().equals(transcriptToFind))
        {
            return variant.canonicalImpact();
        }

        for(PurpleTranscriptImpact otherImpact : variant.otherImpacts())
        {
            if(otherImpact.transcript().equals(transcriptToFind))
            {
                return otherImpact;
            }
        }

        return null;
    }

    @NotNull
    @VisibleForTesting
    static String determineImpact(@NotNull PurpleTranscriptImpact impact)
    {
        String hgvsProteinImpact = impact.hgvsProteinImpact();
        if(!hgvsProteinImpact.isEmpty() && !hgvsProteinImpact.equals("p.?"))
        {
            return AminoAcids.forceSingleLetterProteinAnnotation(hgvsProteinImpact);
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
            joiner.add(VariantEffect.valueOf(effect.name()).effect());
        }
        return joiner.toString();
    }
}
