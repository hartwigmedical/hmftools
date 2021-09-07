package com.hartwig.hmftools.serve.sources.vicc;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableAssociation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotype;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;
import com.hartwig.hmftools.vicc.datamodel.civic.Civic;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivic;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicCoordinates;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicDisease;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicEvidenceItem;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicLifecycleActions;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicProvisionalValue;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicPublicationDate;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicSource;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccTestFactory {

    private static final ViccSource TEST_SOURCE = ViccSource.CGI;

    private ViccTestFactory() {
    }

    @NotNull
    public static ViccEntry testEntryWithGeneEventAndAssociation(@Nullable String gene, @NotNull String event,
            @NotNull Association association) {
        List<Feature> features = Lists.newArrayList(testFeatureWithGeneAndName(gene, event));
        return testViccEntry(TEST_SOURCE, null, features, association);
    }

    @NotNull
    public static ViccEntry testEntryWithOncogenic(@NotNull String oncogenic) {
        return testViccEntry(TEST_SOURCE, null, Lists.newArrayList(), testAssociationWithOncogenic(oncogenic));
    }

    @NotNull
    public static ViccEntry testEntryWithSourceAndTranscript(@NotNull ViccSource source, @Nullable String transcriptId) {
        return testViccEntry(source, transcriptId, Lists.newArrayList(), testAssociation());
    }

    @NotNull
    public static ViccEntry testEntryWithCivicEvidenceDirection(@Nullable String direction) {
        Civic civic = ImmutableCivic.builder()
                .entrezId(Strings.EMPTY)
                .entrezName(Strings.EMPTY)
                .name(Strings.EMPTY)
                .type(Strings.EMPTY)
                .coordinates(ImmutableCivicCoordinates.builder().build())
                .evidenceItem(ImmutableCivicEvidenceItem.builder()
                        .name(Strings.EMPTY)
                        .type(Strings.EMPTY)
                        .status(Strings.EMPTY)
                        .rating(Strings.EMPTY)
                        .evidenceType(Strings.EMPTY)
                        .evidenceLevel(Strings.EMPTY)
                        .evidenceDirection(direction)
                        .disease(ImmutableCivicDisease.builder()
                                .name(Strings.EMPTY)
                                .displayName(Strings.EMPTY)
                                .url(Strings.EMPTY)
                                .id(Strings.EMPTY)
                                .build())
                        .source(ImmutableCivicSource.builder()
                                .status(Strings.EMPTY)
                                .citation(Strings.EMPTY)
                                .sourceUrl(Strings.EMPTY)
                                .pubmedId(Strings.EMPTY)
                                .isReview(Strings.EMPTY)
                                .publicationDate(ImmutableCivicPublicationDate.builder().build())
                                .id(Strings.EMPTY)
                                .build())
                        .openChangeCount(Strings.EMPTY)
                        .description(Strings.EMPTY)
                        .id(Strings.EMPTY)
                        .build())
                .civicActionabilityScore(Strings.EMPTY)
                .alleleRegistryId(Strings.EMPTY)
                .provisionalValue(ImmutableCivicProvisionalValue.builder().build())
                .lifecycleActions(ImmutableCivicLifecycleActions.builder().build())
                .id(Strings.EMPTY)
                .geneId(Strings.EMPTY)
                .description(Strings.EMPTY)
                .build();

        return testViccEntry(ViccSource.CIVIC, null, Lists.newArrayList(), testAssociation(), civic);
    }

    @NotNull
    public static ViccEntry testViccEntry(@NotNull ViccSource source, @Nullable String transcriptId, @NotNull List<Feature> features,
            @NotNull Association association) {
        return testViccEntry(source, transcriptId, features, association, new TestKbSpecificObject());
    }

    @NotNull
    public static ViccEntry testViccEntry(@NotNull ViccSource source, @Nullable String transcriptId, @NotNull List<Feature> features,
            @NotNull Association association, @NotNull KbSpecificObject kbSpecificObject) {
        return ImmutableViccEntry.builder()
                .source(source)
                .association(association)
                .transcriptId(transcriptId)
                .kbSpecificObject(kbSpecificObject)
                .features(features)
                .build();
    }

    @NotNull
    public static Feature testFeatureWithName(@NotNull String name) {
        return testFeatureWithGeneAndName("any", name);
    }

    @NotNull
    public static Feature testFeatureWithGeneAndName(@Nullable String geneSymbol, @NotNull String name) {
        return ImmutableFeature.builder().geneSymbol(geneSymbol).name(name).build();
    }

    @NotNull
    public static Association testActionableAssociation(@Nullable String drugLabels, @NotNull String phenotypeDescription,
            @NotNull String phenotypeTypeId, @Nullable String evidenceLabel, @Nullable String responseType, @Nullable String publication) {
        EvidenceInfo evidenceInfo = null;
        if (publication != null) {
            evidenceInfo = ImmutableEvidenceInfo.builder().addPublications(publication).build();
        }

        EvidenceType evidenceType = ImmutableEvidenceType.builder().sourceName(Strings.EMPTY).build();
        Evidence evidence = ImmutableEvidence.builder().evidenceType(evidenceType).info(evidenceInfo).build();

        Phenotype phenotype = ImmutablePhenotype.builder()
                .description(phenotypeDescription)
                .family(Strings.EMPTY)
                .type(ImmutablePhenotypeType.builder().source(Strings.EMPTY).term(Strings.EMPTY).id(phenotypeTypeId).build())
                .build();

        return ImmutableAssociation.builder()
                .evidence(evidence)
                .evidenceLevel(Strings.EMPTY)
                .evidenceLabel(evidenceLabel)
                .responseType(responseType)
                .drugLabels(drugLabels)
                .sourceLink(Strings.EMPTY)
                .phenotype(phenotype)
                .description(Strings.EMPTY)
                .oncogenic(Strings.EMPTY)
                .build();
    }

    @NotNull
    private static Association testAssociation() {
        return testAssociationWithOncogenic(Strings.EMPTY);
    }

    @NotNull
    private static Association testAssociationWithOncogenic(@NotNull String oncogenic) {
        EvidenceType evidenceType = ImmutableEvidenceType.builder().sourceName(Strings.EMPTY).build();
        Evidence evidence = ImmutableEvidence.builder().evidenceType(evidenceType).build();

        Phenotype phenotype = ImmutablePhenotype.builder().description(Strings.EMPTY).family(Strings.EMPTY).build();

        return ImmutableAssociation.builder()
                .evidence(evidence)
                .evidenceLevel(Strings.EMPTY)
                .evidenceLabel(Strings.EMPTY)
                .responseType(Strings.EMPTY)
                .drugLabels(Strings.EMPTY)
                .sourceLink(Strings.EMPTY)
                .phenotype(phenotype)
                .description(Strings.EMPTY)
                .oncogenic(oncogenic)
                .build();
    }

    private static class TestKbSpecificObject implements KbSpecificObject {

    }
}
