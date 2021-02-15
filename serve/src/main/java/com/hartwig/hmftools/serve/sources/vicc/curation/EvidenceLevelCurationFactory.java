package com.hartwig.hmftools.serve.sources.vicc.curation;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;

final class EvidenceLevelCurationFactory {

    static final Map<EvidenceLevelCurationKey, EvidenceLevel> EVIDENCE_LEVEL_MAPPINGS = Maps.newHashMap();

    private EvidenceLevelCurationFactory() {
    }

    static {
        // This evidence is plain wrong
        EVIDENCE_LEVEL_MAPPINGS.put(trastuzumabBLevelPIK3CA(), EvidenceLevel.D);

        // This evidence is defined pan-cancer which is incorrect
        EVIDENCE_LEVEL_MAPPINGS.put(olaparibALevelBRCA1(), EvidenceLevel.B);
        EVIDENCE_LEVEL_MAPPINGS.put(olaparibALevelBRCA2(), EvidenceLevel.B);

        // This key blocks several evidences of which at least some seem dubious
        EVIDENCE_LEVEL_MAPPINGS.put(everolimusBLevelPTEN(), EvidenceLevel.C);
    }

    @NotNull
    private static EvidenceLevelCurationKey trastuzumabBLevelPIK3CA() {
        return ImmutableEvidenceLevelCurationKey.builder()
                .source(ViccSource.CIVIC)
                .genes(Lists.newArrayList("PIK3CA"))
                .treatment("Trastuzumab")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .build();
    }

    @NotNull
    private static EvidenceLevelCurationKey olaparibALevelBRCA1() {
        return ImmutableEvidenceLevelCurationKey.builder()
                .source(ViccSource.CIVIC)
                .genes(Lists.newArrayList("BRCA1"))
                .treatment("Olaparib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .build();
    }

    @NotNull
    private static EvidenceLevelCurationKey olaparibALevelBRCA2() {
        return ImmutableEvidenceLevelCurationKey.builder()
                .source(ViccSource.CIVIC)
                .genes(Lists.newArrayList("BRCA2"))
                .treatment("Olaparib")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .build();
    }

    @NotNull
    private static EvidenceLevelCurationKey everolimusBLevelPTEN() {
        return ImmutableEvidenceLevelCurationKey.builder()
                .source(ViccSource.CIVIC)
                .genes(Lists.newArrayList("PTEN"))
                .treatment("Everolimus")
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .build();
    }
}
