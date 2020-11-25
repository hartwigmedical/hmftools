package com.hartwig.hmftools.serve.sources.vicc;

import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableAssociation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotype;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKb;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbGene;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbVariant;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKb;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbGene;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbVariant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccTestFactory {

    private ViccTestFactory() {
    }

    @NotNull
    public static ViccEntry testViccEntryWithOncogenic(@NotNull String oncogenic) {
        return testViccEntry(ViccSource.ONCOKB, oncogenic, null);
    }

    @NotNull
    public static ViccEntry testViccEntryWithSourceAndKbObject(@NotNull ViccSource source, @Nullable String transcriptId) {
        return testViccEntry(source, Strings.EMPTY, transcriptId);
    }

    @NotNull
    public static ViccEntry testViccEntry(@NotNull ViccSource source, @NotNull String oncogenic, @Nullable String transcriptId) {
        EvidenceType evidenceType = ImmutableEvidenceType.builder().sourceName(Strings.EMPTY).build();
        Evidence evidence = ImmutableEvidence.builder().evidenceType(evidenceType).build();

        Phenotype phenotype = ImmutablePhenotype.builder().description(Strings.EMPTY).family(Strings.EMPTY).build();

        Association association = ImmutableAssociation.builder()
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

        return ImmutableViccEntry.builder()
                .source(source)
                .association(association)
                .transcriptId(transcriptId)
                .kbSpecificObject(testOncoKb())
                .build();
    }

    @NotNull
    public static Feature testFeatureWithName(@NotNull String name) {
        return testFeatureWithGeneAndName(null, name);
    }

    @NotNull
    public static Feature testFeatureWithGeneAndName(@Nullable String geneSymbol, @NotNull String name) {
        return ImmutableFeature.builder().geneSymbol(geneSymbol).name(name).build();
    }

    @NotNull
    private static OncoKb testOncoKb() {
        OncoKbConsequence consequence = ImmutableOncoKbConsequence.builder()
                .term(Strings.EMPTY)
                .description(Strings.EMPTY)
                .isGenerallyTruncating(Strings.EMPTY)
                .build();

        OncoKbGene gene = ImmutableOncoKbGene.builder()
                .hugoSymbol(Strings.EMPTY)
                .name(Strings.EMPTY)
                .entrezGeneId(Strings.EMPTY)
                .curatedIsoform(Strings.EMPTY)
                .curatedRefSeq(Strings.EMPTY)
                .oncogene(Strings.EMPTY)
                .tsg(Strings.EMPTY)
                .build();

        OncoKbVariant variant = ImmutableOncoKbVariant.builder()
                .name(Strings.EMPTY)
                .alteration(Strings.EMPTY)
                .consequence(consequence)
                .gene(gene)
                .proteinStart(Strings.EMPTY)
                .proteinEnd(Strings.EMPTY)
                .refResidues(Strings.EMPTY)
                .variantResidues(Strings.EMPTY)
                .build();

        OncoKbBiological biological = ImmutableOncoKbBiological.builder()
                .gene(Strings.EMPTY)
                .entrezGeneId(Strings.EMPTY)
                .isoform(Strings.EMPTY)
                .refSeq(Strings.EMPTY)
                .oncokbVariant(variant)
                .oncogenic(Strings.EMPTY)
                .mutationEffect(Strings.EMPTY)
                .mutationEffectPmids(Strings.EMPTY)
                .mutationEffectAbstracts(Strings.EMPTY)
                .build();

        return ImmutableOncoKb.builder().oncoKbBiological(biological).build();
    }
}
