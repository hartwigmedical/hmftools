package com.hartwig.hmftools.ckb.interpretation.common;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodelinterpretation.common.ImmutableReferenceExtend;
import com.hartwig.hmftools.ckb.datamodelinterpretation.common.ReferenceExtend;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.reference.Reference;

import org.jetbrains.annotations.NotNull;

public class CommonInterpretationFactory {

    private CommonInterpretationFactory() {

    }

    @NotNull
    public static List<ReferenceExtend> extractReferences(@NotNull List<ReferenceInfo> referenceInfos, @NotNull CkbJsonDatabase ckbEntry) {
        List<ReferenceExtend> references = Lists.newArrayList();
        for (ReferenceInfo referenceInfo : referenceInfos) {
            for (Reference reference : ckbEntry.references()) {
                if (referenceInfo.id() == reference.id()) {
                    references.add(ImmutableReferenceExtend.builder()
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
}
