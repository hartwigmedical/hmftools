package com.hartwig.hmftools.orange.report.datamodel;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SmallVariantFactory
{
    @NotNull
    public static List<SmallVariant> create(@NotNull List<PurpleVariant> variants, @NotNull List<PurpleDriver> drivers)
    {
        List<SmallVariant> entries = new ArrayList<>();
        for(PurpleVariant variant : variants)
        {
            if(hasCanonicalImpact(variant))
            {
                PurpleDriver driver = Drivers.canonicalMutationEntryForGene(drivers, variant.gene());
                entries.add(toVariantEntry(variant, driver));
            }
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
    private static SmallVariant toVariantEntry(@NotNull PurpleVariant variant, @Nullable PurpleDriver driver)
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

        return ImmutableSmallVariant.builder()
                .purpleVariant(variant)
                .driver(driver)
                .transcriptImpact(transcriptImpact)
                .isCanonical(driver == null || driver.transcript().equals(variant.canonicalImpact().transcript()))
                .driverInterpretation(driver != null ? DriverInterpretation.interpret(driver.driverLikelihood()) : null)
                .build();
    }

    private static String findingKey()
    {

    }

    @NotNull
    private static List<PurpleVariant> findReportedVariantsForDriver(@NotNull List<PurpleVariant> variants, @NotNull PurpleDriver driver)
    {
        List<PurpleVariant> reportedVariantsForDriver = new ArrayList<>();
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
        List<PurpleVariant> reportedVariantsForGene = new ArrayList<>();
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

    private static boolean hasCanonicalImpact(@NotNull PurpleVariant variant)
    {
        return !variant.canonicalImpact().transcript().isEmpty();
    }
}
