package com.hartwig.hmftools.ckb.datamodel.common;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.indication.JsonIndication;
import com.hartwig.hmftools.ckb.json.reference.JsonReference;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CommonInterpretationFactory {

    private CommonInterpretationFactory() {
    }

    @NotNull
    public static List<com.hartwig.hmftools.ckb.datamodel.common.Reference> extractReferences(
            @NotNull List<ReferenceInfo> referenceInfos, @NotNull CkbJsonDatabase ckbEntry) {
        List<com.hartwig.hmftools.ckb.datamodel.common.Reference> references = Lists.newArrayList();
        for (ReferenceInfo referenceInfo : referenceInfos) {
            for (JsonReference reference : ckbEntry.references()) {
                if (referenceInfo.id() == reference.id()) {
                    references.add(ImmutableReference.builder()
                            .id(reference.id())
                            .pubMedId(reference.pubMedId())
                            .title(reference.title())
                            .url(reference.url())
                            .authors(reference.authors())
                            .journal(reference.journal())
                            .volume(reference.volume())
                            .issue(reference.issue())
                            .date(reference.date())
                            .abstractText(reference.abstractText())
                            .year(reference.year())
                            .build());
                }
            }
        }
        return references;
    }

    @NotNull
    public static com.hartwig.hmftools.ckb.datamodel.common.Indication extractIndication(
            @NotNull CkbJsonDatabase ckbEntry, @Nullable IndicationInfo indicationInfo) {
        ImmutableIndication.Builder outputBuilder = ImmutableIndication.builder();
        for (JsonIndication indication : ckbEntry.indications()) {
            if (indicationInfo.id().equals(indication.id())) {
                outputBuilder.id(indication.id())
                        .name(indication.name())
                        .source(indication.source())
                        .definition(indication.definition())
                        .currentPreferredTerm(indication.currentPreferredTerm())
                        .lastUpdateDateFromDO(indication.lastUpdateDateFromDO())
                        .altIds(indication.altIds())
                        .termId(indication.termId());
            }
        }
        return outputBuilder.build();
    }
}
