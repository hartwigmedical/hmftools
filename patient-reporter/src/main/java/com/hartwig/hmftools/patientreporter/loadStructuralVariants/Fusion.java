package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class Fusion {

    public abstract boolean reportable();

    @NotNull
    public abstract String knownType();

    @NotNull
    public abstract String primarySource();

    @NotNull
    public abstract String clusterId();

    @NotNull
    public abstract String clusterCount();

    @NotNull
    public abstract String resolvedType();

    @NotNull
    public abstract String svIdUp();

    @NotNull
    public abstract String chrUp();

    @NotNull
    public abstract String posUp();

    @NotNull
    public abstract String orientUp();

    @NotNull
    public abstract String typeUp();

    public abstract double ploidyUp();

    @NotNull
    public abstract String geneUp();

    @NotNull
    public abstract String chrBandUp();

    @NotNull
    public abstract String transcriptUp();

    @NotNull
    public abstract String strandUp();

    @NotNull
    public abstract String regionTypeUp();

    @NotNull
    public abstract String codingTypeUp();

    public abstract int exonUp();

    @NotNull
    public abstract String phaseUp();

    @NotNull
    public abstract String exonMaxUp();

    @NotNull
    public abstract String disruptiveUp();

    @NotNull
    public abstract String exactBaseUp();

    @NotNull
    public abstract String codingBasesUp();

    @NotNull
    public abstract String totalCodingUp();

    @NotNull
    public abstract String codingStartUp();

    @NotNull
    public abstract String codingEndUp();

    @NotNull
    public abstract String transStartUp();

    @NotNull
    public abstract String transEndUp();

    @NotNull
    public abstract String distancePrevUp();

    @NotNull
    public abstract String canonicalUp();

    @NotNull
    public abstract String biotypeUp();

    @NotNull
    public abstract String svIdDown();

    @NotNull
    public abstract String chrDown();

    @NotNull
    public abstract String posDown();

    @NotNull
    public abstract String orientDown();

    @NotNull
    public abstract String typeDown();

    public abstract double ploidyDown();

    @NotNull
    public abstract String geneDown();

    @NotNull
    public abstract String chrBandDown();

    @NotNull
    public abstract String transcriptDown();

    @NotNull
    public abstract String strandDown();

    @NotNull
    public abstract String regionTypeDown();

    @NotNull
    public abstract String codingTypeDown();

    public abstract int exonDown();

    @NotNull
    public abstract String phaseDown();

    @NotNull
    public abstract String exonMaxDown();

    @NotNull
    public abstract String disruptiveDown();

    @NotNull
    public abstract String exactBaseDown();

    @NotNull
    public abstract String codingBasesDown();

    @NotNull
    public abstract String totalCodingDown();

    @NotNull
    public abstract String codingStartDown();

    @NotNull
    public abstract String codingEndDown();

    @NotNull
    public abstract String transStartDown();

    @NotNull
    public abstract String transEndDown();

    @NotNull
    public abstract String distancePrevDown();

    @NotNull
    public abstract String canonicalDown();

    @NotNull
    public abstract String biotypeDown();

    @NotNull
    public abstract String proteinsKept();

    @NotNull
    public abstract String proteinsLost();
}
