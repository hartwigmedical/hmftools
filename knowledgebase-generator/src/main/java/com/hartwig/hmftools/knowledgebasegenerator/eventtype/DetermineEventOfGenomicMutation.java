package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.knowledgebasegenerator.GenomicEvents;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.ActionableAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.CnvExtractor;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.ImmutableActionableAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.ImmutableKnownAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.KnownAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.fusion.ImmutableKnownFusions;
import com.hartwig.hmftools.knowledgebasegenerator.fusion.KnownFusions;
import com.hartwig.hmftools.knowledgebasegenerator.signatures.ImmutableSignatures;
import com.hartwig.hmftools.knowledgebasegenerator.signatures.Signatures;
import com.hartwig.hmftools.knowledgebasegenerator.signatures.SignaturesExtractor;
import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class DetermineEventOfGenomicMutation {
    private static final Logger LOGGER = LogManager.getLogger(DetermineEventOfGenomicMutation.class);

    private static final List<String> AMPLIFICATION = Lists.newArrayList("Amplification", "AMPLIFICATION", "amplification");

    private static final List<String> DELETION = Lists.newArrayList("Deletion", "DELETION", "deletion");

    private static final Set<String> VARIANTS = Sets.newHashSet("missense_variant", "inframe_deletion", "inframe_insertion");
    private static final Set<String> FUSIONS = Sets.newHashSet("Fusion", "Fusions", "FUSIONS", "Gene Fusion", "Transcript Fusion");
    private static final Set<String> RANGE = Sets.newHashSet();
    private static final Set<String> SIGNATURE_MSI = Sets.newHashSet("Microsatellite Instability-High");
    private static final Set<String> SIGNATURE_HRD = Sets.newHashSet("");
    private static final Set<String> SIGNATURE_MTL = Sets.newHashSet("");
    private static final Set<String> SIGNATURE_MTB = Sets.newHashSet("");

    @NotNull
    public static KnownAmplificationDeletion checkKnownAmplification(@NotNull ViccEntry viccEntry, @NotNull EventType type,
            @NotNull String gene, @NotNull String event) {
        Source source = Source.sourceFromKnowledgebase(viccEntry.source());
        if (AMPLIFICATION.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Amplification");
            return CnvExtractor.determineKnownAmplificationDeletion(source, typeEvent.toString(), gene);
        }
        return ImmutableKnownAmplificationDeletion.builder().gene("").eventType("").source("").sourceLink("").build();
    }

    @NotNull
    public static KnownAmplificationDeletion checkKnownDeletion(@NotNull ViccEntry viccEntry, @NotNull EventType type, @NotNull String gene,
            @NotNull String event) {
        Source source = Source.sourceFromKnowledgebase(viccEntry.source());

        if (DELETION.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Deletion");
            return CnvExtractor.determineKnownAmplificationDeletion(source, typeEvent.toString(), gene);
        }
        return ImmutableKnownAmplificationDeletion.builder().gene("").eventType("").source("").sourceLink("").build();
    }

    @NotNull
    public static ActionableAmplificationDeletion checkActionableAmplification(@NotNull ViccEntry viccEntry, @NotNull EventType type,
            @NotNull String gene, @NotNull String event) {
        Source source = Source.sourceFromKnowledgebase(viccEntry.source());
        if (AMPLIFICATION.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Amplification");
            return CnvExtractor.determineActionableAmplificationDeletion(source, typeEvent.toString(), gene, viccEntry);

        }
        return ImmutableActionableAmplificationDeletion.builder()
                .gene("")
                .eventType("")
                .source("")
                .drug("")
                .drugType("")
                .cancerType("")
                .level("")
                .direction("")
                .sourceLink("")
                .build();
    }

    @NotNull
    public static ActionableAmplificationDeletion checkActionableDeletion(@NotNull ViccEntry viccEntry, @NotNull EventType type,
            @NotNull String gene, @NotNull String event) {
        Source source = Source.sourceFromKnowledgebase(viccEntry.source());

        if (DELETION.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Deletion");
            return CnvExtractor.determineActionableAmplificationDeletion(source, typeEvent.toString(), gene, viccEntry);
        }
        return ImmutableActionableAmplificationDeletion.builder()
                .gene("")
                .eventType("")
                .source("")
                .drug("")
                .drugType("")
                .cancerType("")
                .level("")
                .direction("")
                .sourceLink("")
                .build();
    }

    public static void checkVariants(@NotNull ViccEntry viccEntry, @NotNull EventType type, @NotNull String gene, @NotNull String event) {
        Source source = Source.sourceFromKnowledgebase(viccEntry.source());

        if (VARIANTS.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Variants");
        }
    }

    public static void checkRange(@NotNull ViccEntry viccEntry, @NotNull EventType type, @NotNull String gene, @NotNull String event) {
        Source source = Source.sourceFromKnowledgebase(viccEntry.source());

        if (RANGE.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Range");
        }
    }

    @NotNull
    public static KnownFusions checkFusions(@NotNull ViccEntry viccEntry, @NotNull EventType type, @NotNull String gene,
            @NotNull String event) {
        Source source = Source.sourceFromKnowledgebase(viccEntry.source());

        if (FUSIONS.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Fusions");
        }
        return ImmutableKnownFusions.builder().gene("").eventType("").source("").sourceLink("").build();
    }

    @NotNull
    public static Signatures checkSignatures(@NotNull ViccEntry viccEntry, @NotNull EventType type, @NotNull String gene, @NotNull String event) {
        Source source = Source.sourceFromKnowledgebase(viccEntry.source());
            if (SIGNATURE_MSI.contains(event)){
                GenomicEvents typeEvent = GenomicEvents.genomicEvents("MSI");
                return SignaturesExtractor.determineSignatures(source, typeEvent.toString(), gene);
            } else if (SIGNATURE_HRD.contains(event)) {
                GenomicEvents typeEvent = GenomicEvents.genomicEvents("HRD");
                return SignaturesExtractor.determineSignatures(source, typeEvent.toString(), gene);
            } else if (SIGNATURE_MTL.contains(event)) {
                GenomicEvents typeEvent = GenomicEvents.genomicEvents("MTL");
                return SignaturesExtractor.determineSignatures(source, typeEvent.toString(), gene);
            } else if (SIGNATURE_MTB.contains(event)) {
                GenomicEvents typeEvent = GenomicEvents.genomicEvents("MTB");
                return SignaturesExtractor.determineSignatures(source, typeEvent.toString(), gene);
            }

        return ImmutableSignatures.builder().eventType(Strings.EMPTY).source(Strings.EMPTY).sourceLink(Strings.EMPTY).build();
    }
}
