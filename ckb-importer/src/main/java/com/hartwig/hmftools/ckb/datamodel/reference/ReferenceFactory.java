package com.hartwig.hmftools.ckb.datamodel.reference;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.reference.JsonReference;

import org.jetbrains.annotations.NotNull;

public final class ReferenceFactory {

    private ReferenceFactory() {
    }

    @NotNull
    public static List<Reference> extractReferences(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull List<ReferenceInfo> referenceInfos) {
        List<Reference> references = Lists.newArrayList();
        for (ReferenceInfo referenceInfo : referenceInfos) {
            references.add(resolveReference(ckbJsonDatabase, referenceInfo));
        }
        return references;
    }

    @NotNull
    private static Reference resolveReference(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull ReferenceInfo referenceInfo) {
        for (JsonReference reference : ckbJsonDatabase.references()) {
            if (reference.id() == referenceInfo.id()) {
                return ImmutableReference.builder()
                        .id(reference.id())
                        .pubMedId(reference.pubMedId())
                        .title(reference.title())
                        .abstractText(reference.abstractText())
                        .url(reference.url())
                        .journal(reference.journal())
                        .authors(reference.authors())
                        .volume(reference.volume())
                        .issue(reference.issue())
                        .date(reference.date())
                        .year(reference.year())
                        .build();
            }
        }

        throw new IllegalStateException("Could not resolve CKB reference with id '" + referenceInfo.id() + "'");
    }
}
