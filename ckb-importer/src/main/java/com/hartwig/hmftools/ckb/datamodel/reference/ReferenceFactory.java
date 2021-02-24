package com.hartwig.hmftools.ckb.datamodel.reference;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.reference.JsonReference;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReferenceFactory {

    private static final Logger LOGGER = LogManager.getLogger(ReferenceFactory.class);

    private ReferenceFactory() {
    }

    @Nullable
    public static String extractDescription(@NotNull String type, int id, @NotNull List<DescriptionInfo> descriptionInfos) {
        if (descriptionInfos.isEmpty()) {
            return null;
        }

        if (descriptionInfos.size() > 1) {
            LOGGER.warn("Multiple descriptions found for object of type '{}' with id '{}'", type, id);
        }

        return descriptionInfos.get(0).description();
    }

    @NotNull
    public static List<Reference> extractDescriptionReferences(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<DescriptionInfo> descriptionInfos) {
        List<Reference> references = Lists.newArrayList();
        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            references.addAll(extractReferences(ckbJsonDatabase, descriptionInfo.references()));
        }
        return references;
    }

    @NotNull
    public static List<Reference> extractReferences(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull List<ReferenceInfo> referenceInfos) {
        List<Reference> references = Lists.newArrayList();
        for (ReferenceInfo referenceInfo : referenceInfos) {
            Reference resolvedReference = resolveReference(ckbJsonDatabase, referenceInfo);
            // References that miss both title and pubmed are generic references to certain websites and not informative, so can leave out.
            if (resolvedReference.title() != null || resolvedReference.pubMedId() != null) {
                references.add(resolvedReference);
            }
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
