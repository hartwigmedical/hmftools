package com.hartwig.hmftools.datamodel.finding;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.finding.clinicaltranscript.ClinicalTranscriptsModel;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class SmallVariantFactory
{
    public static DriverFindings<SmallVariant> smallVariantFindings(@NotNull PurpleRecord purpleRecord, @Nullable ClinicalTranscriptsModel clinicalTranscriptsModel,
            @NotNull Map<String, DriverGene> driverGeneMap) {
        List<SmallVariant> allSmallVariants = Lists.newArrayList();
        allSmallVariants.addAll(SmallVariantFactory.create(
                DriverSource.SOMATIC, purpleRecord.reportableSomaticVariants(), purpleRecord.somaticDrivers(),
                clinicalTranscriptsModel, driverGeneMap));

        List<PurpleVariant> germlineVariants = purpleRecord.reportableGermlineVariants();
        List<PurpleDriver> germlineDrivers = purpleRecord.germlineDrivers();
        if (germlineVariants != null && germlineDrivers != null)
        {
            allSmallVariants.addAll(SmallVariantFactory.create(
                    DriverSource.GERMLINE, germlineVariants, germlineDrivers, clinicalTranscriptsModel, driverGeneMap));
        }
        return ImmutableDriverFindings.<SmallVariant>builder()
                .findingsStatus(purpleStatus(purpleRecord))
                .findings(allSmallVariants)
                .build();
    }

    private static FindingsStatus purpleStatus(PurpleRecord purpleRecord) {
        return purpleRecord.fit().qc().status().equals(Set.of(PurpleQCStatus.PASS)) ? FindingsStatus.OK : FindingsStatus.NOT_AVAILABLE;
    }

    @NotNull
    private static List<SmallVariant> create(
            @NotNull DriverSource sampleType, @NotNull List<PurpleVariant> variants, @NotNull List<PurpleDriver> drivers,
            @Nullable ClinicalTranscriptsModel clinicalTranscriptsModel,
            @NotNull Map<String, DriverGene> driverGeneMap)
    {
        List<SmallVariant> entries = new ArrayList<>();
        for(PurpleVariant variant : variants)
        {
            if(hasCanonicalImpact(variant))
            {
                // find the driver object, if it is found, it is reported finding
                PurpleDriver driver = Drivers.canonicalMutationEntryForGene(drivers, variant.gene());
                if(driver != null)
                {
                    entries.add(toSmallVariant(variant, driver, sampleType, clinicalTranscriptsModel, driverGeneMap));
                }
            }
        }

        for(PurpleDriver nonCanonicalDriver : Drivers.nonCanonicalMutationEntries(drivers))
        {
            List<PurpleVariant> nonCanonicalVariants = findReportedVariantsForDriver(variants, nonCanonicalDriver);
            for(PurpleVariant nonCanonicalVariant : nonCanonicalVariants)
            {
                entries.add(toSmallVariant(nonCanonicalVariant, nonCanonicalDriver, sampleType, clinicalTranscriptsModel, driverGeneMap));
            }
        }

        return entries;
    }

    @NotNull
    private static SmallVariant toSmallVariant(@NotNull PurpleVariant variant, @NotNull PurpleDriver driver,
            @NotNull DriverSource sampleType, @Nullable ClinicalTranscriptsModel clinicalTranscriptsModel, @NotNull Map<String, DriverGene> driverGeneMap)
    {
        PurpleTranscriptImpact transcriptImpact;

        transcriptImpact = findTranscriptImpact(variant, driver.transcript());
        if(transcriptImpact == null)
        {
            throw new IllegalStateException("Could not find impact on transcript " + driver.transcript() + " for variant " + variant);
        }

        boolean isCanonical = driver.transcript().equals(variant.canonicalImpact().transcript());

        PurpleTranscriptImpact otherImpact = isCanonical && clinicalTranscriptsModel != null ? findOtherImpactClinical(variant, clinicalTranscriptsModel)
                : null;

        DriverGene driverGene = driverGeneMap.get(variant.gene());
        DriverCategory driverCategory = driverGene != null ? driverLikelihoodType(driverGene.likelihoodType()) : null;
        return ImmutableSmallVariant.builder()
                .findingKey(FindingKeys.smallVariant(sampleType, variant, transcriptImpact, isCanonical))
                .driverSource(sampleType)
                .driverLikelihoodType(driverCategory)
                .reportedStatus(ReportedStatus.REPORTED) // all drivers here are reported
                .driverInterpretation(DriverInterpretation.interpret(driver.driverLikelihood()))
                .purpleVariant(variant)
                .driver(driver)
                .transcriptImpact(transcriptImpact)
                .otherImpact(otherImpact)
                .isCanonical(isCanonical)
                .build();
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

        return findOtherTranscriptImpact(variant, transcriptToFind);
    }

    @Nullable
    static PurpleTranscriptImpact findOtherTranscriptImpact(@NotNull PurpleVariant variant, @NotNull String transcriptToFind)
    {
        // ADQ: Recommended implementation:
        //return variant.otherImpacts().stream().filter(o -> o.transcript().equals(transcriptToFind)).findFirst().orElse(null);

        for(PurpleTranscriptImpact otherImpact : variant.otherImpacts())
        {
            if(otherImpact.transcript().equals(transcriptToFind))
            {
                return otherImpact;
            }
        }

        return null;
    }

    @Nullable
    private static PurpleTranscriptImpact findOtherImpactClinical(@NotNull PurpleVariant variant, @NotNull ClinicalTranscriptsModel clinicalTranscriptsModel)
    {
        String transcriptOverride = clinicalTranscriptsModel.findCanonicalTranscriptForGene(variant.gene());
        if (transcriptOverride != null)
        {
            PurpleTranscriptImpact otherImpactClinical = findOtherTranscriptImpact(variant, transcriptOverride);
            if (otherImpactClinical != null)
            {
                otherImpactClinical = !otherImpactClinical.hgvsCodingImpact().equals(variant.canonicalImpact().hgvsCodingImpact())
                        ? otherImpactClinical
                        : null;
            }
            return otherImpactClinical;
        }
        return null;
    }

    private static boolean hasCanonicalImpact(@NotNull PurpleVariant variant)
    {
        return !variant.canonicalImpact().transcript().isEmpty();
    }

    private static DriverCategory driverLikelihoodType(@NotNull com.hartwig.hmftools.common.driver.DriverCategory driverLikelihoodType) {
        return switch (driverLikelihoodType) {
            case ONCO -> DriverCategory.ONCO;
            case TSG -> DriverCategory.TSG;
        };
    }
}
