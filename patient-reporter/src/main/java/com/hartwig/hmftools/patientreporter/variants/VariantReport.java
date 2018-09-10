package com.hartwig.hmftools.patientreporter.variants;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VariantReport {

    @NotNull
    SomaticVariant variant();

    @NotNull
    String gene();

    @NotNull
    String variantDetails();

    int totalReadCount();

    int alleleReadCount();

    @Value.Derived
    default double alleleFrequency() {
        return (double) alleleReadCount() / totalReadCount();
    }

    @NotNull
    @Value.Derived
    default String readDepth() {
        return alleleReadCount() + " / " + totalReadCount() + " (" + PatientReportFormat.formatPercent(alleleFrequency()) + ")";
    }

    @Nullable
    String cosmicID();

    @NotNull
    @Value.Derived
    default String cosmicUrl() {
        String COSMIC_IDENTIFIER = "COSM";
        String cosmicID = cosmicID();

        if (cosmicID == null) {
            return Strings.EMPTY;
        }
        final int identifierPos = cosmicID.indexOf(COSMIC_IDENTIFIER);
        String cosmicIdentifier;
        if (identifierPos >= 0) {
            cosmicIdentifier = cosmicID.substring(identifierPos + COSMIC_IDENTIFIER.length());
        } else {
            cosmicIdentifier = cosmicID;
        }

        return "http://cancer.sanger.ac.uk/cosmic/mutation/overview?genome=37&id=" + cosmicIdentifier;
    }

    @NotNull
    String ploidy();

    double purityAdjustedVAF();

    @NotNull
    @Value.Derived
    default String ploidyVaf() {
        return ploidy() + " (" + PatientReportFormat.formatPercent(purityAdjustedVAF()) + ")";
    }

    double clonalProbability();

    @NotNull
    String wildTypeStatus();

    double driverProbability();

    @NotNull
    String actionabilityLevel();
}
