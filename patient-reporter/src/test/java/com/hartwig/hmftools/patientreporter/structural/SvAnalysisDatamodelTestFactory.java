package com.hartwig.hmftools.patientreporter.structural;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class SvAnalysisDatamodelTestFactory {

    private SvAnalysisDatamodelTestFactory() {
    }

    @NotNull
    static ImmutableDisruption.Builder disruptionBuilder() {
        return ImmutableDisruption.builder()
                .reportable(false)
                .svId(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .position(Strings.EMPTY)
                .orientation(1)
                .type(Strings.EMPTY)
                .ploidy(1)
                .gene(Strings.EMPTY)
                .chrBand(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .strand(1)
                .regionType(Strings.EMPTY)
                .codingType(Strings.EMPTY)
                .canonical(false)
                .biotype(Strings.EMPTY)
                .exonUp(0)
                .exonDown(0)
                .isDisruptive(false);
    }

    @NotNull
    static ImmutableFusion.Builder fusionBuilder() {
        return ImmutableFusion.builder()
                .reportable(true)
                .knownType(Strings.EMPTY)
                .primarySource(Strings.EMPTY)
                .clusterId(Strings.EMPTY)
                .clusterCount(Strings.EMPTY)
                .resolvedType(Strings.EMPTY)
                .svIdUp(Strings.EMPTY)
                .chrUp(Strings.EMPTY)
                .posUp(Strings.EMPTY)
                .orientUp(Strings.EMPTY)
                .typeUp(Strings.EMPTY)
                .ploidyUp(1D)
                .geneUp(Strings.EMPTY)
                .chrBandUp(Strings.EMPTY)
                .transcriptUp(Strings.EMPTY)
                .strandUp(Strings.EMPTY)
                .regionTypeUp(Strings.EMPTY)
                .codingTypeUp(Strings.EMPTY)
                .exonUp(1)
                .phaseUp(Strings.EMPTY)
                .exonMaxUp(Strings.EMPTY)
                .disruptiveUp(Strings.EMPTY)
                .exactBaseUp(Strings.EMPTY)
                .codingBasesUp(Strings.EMPTY)
                .totalCodingUp(Strings.EMPTY)
                .codingStartUp(Strings.EMPTY)
                .codingEndUp(Strings.EMPTY)
                .transStartUp(Strings.EMPTY)
                .transEndUp(Strings.EMPTY)
                .distancePrevUp(Strings.EMPTY)
                .canonicalUp(Strings.EMPTY)
                .biotypeUp(Strings.EMPTY)
                .svIdDown(Strings.EMPTY)
                .chrDown(Strings.EMPTY)
                .posDown(Strings.EMPTY)
                .orientDown(Strings.EMPTY)
                .typeDown(Strings.EMPTY)
                .ploidyDown(1D)
                .geneDown(Strings.EMPTY)
                .chrBandDown(Strings.EMPTY)
                .transcriptDown(Strings.EMPTY)
                .strandDown(Strings.EMPTY)
                .regionTypeDown(Strings.EMPTY)
                .codingTypeDown(Strings.EMPTY)
                .exonDown(14)
                .phaseDown(Strings.EMPTY)
                .exonMaxDown(Strings.EMPTY)
                .disruptiveDown(Strings.EMPTY)
                .exactBaseDown(Strings.EMPTY)
                .codingBasesDown(Strings.EMPTY)
                .totalCodingDown(Strings.EMPTY)
                .codingStartDown(Strings.EMPTY)
                .codingEndDown(Strings.EMPTY)
                .transStartDown(Strings.EMPTY)
                .transEndDown(Strings.EMPTY)
                .distancePrevDown(Strings.EMPTY)
                .canonicalDown(Strings.EMPTY)
                .biotypeDown(Strings.EMPTY)
                .proteinsKept(Strings.EMPTY)
                .proteinsLost(Strings.EMPTY);
    }
}
