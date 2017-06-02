package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.formatPercent;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VariantReport extends GenomePosition {

    @NotNull
    String gene();

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    String transcript();

    @NotNull
    String hgvsCoding();

    @NotNull
    String hgvsProtein();

    @NotNull
    String consequence();

    @NotNull
    @Value.Default
    default String cosmicID() {
        return Strings.EMPTY;
    }

    @NotNull
    @Value.Default
    default String baf() {
        return Strings.EMPTY;
    }

    @NotNull
    @Value.Default
    default double impliedVAF() {
        return 0;
    }

    default String variantAlleleFrequency() {
        return Doubles.isZero(impliedVAF())
                ? formatPercent(alleleFrequency())
                : formatPercent(alleleFrequency()) + " (" + formatPercent(impliedVAF()) + ")";
    }

    int totalReadCount();

    int alleleReadCount();

    default double alleleFrequency() {
        return (double) alleleReadCount() / totalReadCount();
    }
}
