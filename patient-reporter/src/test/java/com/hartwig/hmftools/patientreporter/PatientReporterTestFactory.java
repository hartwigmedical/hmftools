package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.variants.germline.ImmutableGermlineVariant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PatientReporterTestFactory {

    private PatientReporterTestFactory() {
    }

    @NotNull
    public static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder() {
        return ImmutableGeneCopyNumber.builder()
                .start(1)
                .end(2)
                .gene(Strings.EMPTY)
                .chromosome("1")
                .chromosomeBand(Strings.EMPTY)
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .somaticRegions(1)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .transcriptID("trans")
                .transcriptVersion(0)
                .minMinorAllelePloidy(0);
    }

    @NotNull
    public static ImmutableEnrichedSomaticVariant.Builder createTestEnrichedSomaticVariantBuilder() {
        return ImmutableEnrichedSomaticVariant.builder()
                .trinucleotideContext(Strings.EMPTY)
                .highConfidenceRegion(false)
                .microhomology(Strings.EMPTY)
                .repeatSequence(Strings.EMPTY)
                .repeatCount(0)
                .clonality(Clonality.UNKNOWN)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .type(VariantType.UNDEFINED)
                .filter("PASS")
                .totalReadCount(0)
                .alleleReadCount(0)
                .gene(Strings.EMPTY)
                .genesEffected(0)
                .worstEffect(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONE)
                .worstEffectTranscript(Strings.EMPTY)
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .hotspot(Hotspot.NON_HOTSPOT)
                .recovered(false)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAllelePloidy(0)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .ploidy(0)
                .mappability(0);
    }

    @NotNull
    public static ImmutableGermlineVariant.Builder createTestGermlineVariantBuilder() {
        return ImmutableGermlineVariant.builder()
                .passFilter(true)
                .gene(Strings.EMPTY)
                .hgvsCodingImpact(Strings.EMPTY)
                .hgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .germlineStatus(Strings.EMPTY)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAllelePloidy(1D)
                .biallelic(false);
    }
}
