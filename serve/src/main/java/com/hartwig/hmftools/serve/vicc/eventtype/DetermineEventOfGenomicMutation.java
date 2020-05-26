package com.hartwig.hmftools.serve.vicc.eventtype;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.GenomicEvents;
import com.hartwig.hmftools.serve.vicc.copynumber.ActionableAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.copynumber.ImmutableActionableAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.copynumber.ImmutableKnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.fusion.ImmutableKnownFusions;
import com.hartwig.hmftools.serve.vicc.fusion.KnownFusions;
import com.hartwig.hmftools.serve.vicc.signatures.ImmutableSignatures;
import com.hartwig.hmftools.serve.vicc.signatures.Signatures;
import com.hartwig.hmftools.serve.vicc.signatures.SignaturesExtractor;
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
    private static final Set<String> FUSIONS_PAIRS = Sets.newHashSet("fusion pair");
    private static final Set<String> FUSIONS_PROMISCUOUS = Sets.newHashSet("fusion promiscuous");
    private static final Set<String> RANGE = Sets.newHashSet();
    private static final Set<String> SIGNATURE_MSI = Sets.newHashSet("Microsatellite Instability-High");
    private static final Set<String> SIGNATURE_HRD = Sets.newHashSet("");
    private static final Set<String> SIGNATURE_MTL = Sets.newHashSet("");
    private static final Set<String> SIGNATURE_MTB = Sets.newHashSet("");

    @NotNull
    public static KnownAmplificationDeletion checkKnownAmplification(@NotNull ViccEntry viccEntry, @NotNull String gene,
            @NotNull String event) {
        if (AMPLIFICATION.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Amplification");
            return CopyNumberExtractor.determineKnownAmplificationDeletion(viccEntry.source(), typeEvent.toString(), gene);
        }
        return ImmutableKnownAmplificationDeletion.builder().gene("").eventType("").source("").sourceLink("").build();
    }

    @NotNull
    public static KnownAmplificationDeletion checkKnownDeletion(@NotNull ViccEntry viccEntry, @NotNull String gene, @NotNull String event) {
        if (DELETION.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Deletion");
            return CopyNumberExtractor.determineKnownAmplificationDeletion(viccEntry.source(), typeEvent.toString(), gene);
        }
        return ImmutableKnownAmplificationDeletion.builder().gene("").eventType("").source("").sourceLink("").build();
    }

    @NotNull
    public static ActionableAmplificationDeletion checkActionableAmplification(@NotNull ViccEntry viccEntry, @NotNull String gene,
            @NotNull String event) {
        if (AMPLIFICATION.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Amplification");
            return CopyNumberExtractor.determineActionableAmplificationDeletion(viccEntry.source(), typeEvent.toString(), gene, viccEntry);

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
    public static ActionableAmplificationDeletion checkActionableDeletion(@NotNull ViccEntry viccEntry, @NotNull String gene,
            @NotNull String event) {
        if (DELETION.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Deletion");
            return CopyNumberExtractor.determineActionableAmplificationDeletion(viccEntry.source(), typeEvent.toString(), gene, viccEntry);
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

    public static void checkVariants(@NotNull ViccEntry viccEntry, @NotNull String gene, @NotNull String event) {
        if (VARIANTS.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Variants");
        }
    }

    public static void checkRange(@NotNull ViccEntry viccEntry, @NotNull String gene, @NotNull String event) {
        if (RANGE.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("Range");
        }
    }

    @NotNull
    public static KnownFusions checkFusionsPairs(@NotNull ViccEntry viccEntry, @NotNull String gene, @NotNull String event) {
        String function =
                viccEntry.association().evidence().description() == null ? Strings.EMPTY : viccEntry.association().evidence().description();

        if (FUSIONS_PAIRS.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("fusion pair");
            return FusionExtractor.determineKnownFusionsPairs(viccEntry.source(), typeEvent.toString(), gene, function);
        }
        return ImmutableKnownFusions.builder().gene("").eventType("").source("").sourceLink("").build();
    }

    @NotNull
    public static KnownFusions checkFusionPromiscuous(@NotNull ViccEntry viccEntry, @NotNull String gene, @NotNull String event) {
        String function =
                viccEntry.association().evidence().description() == null ? Strings.EMPTY : viccEntry.association().evidence().description();

        if (FUSIONS_PROMISCUOUS.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("fusion promiscuous");
            return FusionExtractor.determinePromiscuousFusions(viccEntry.source(), typeEvent.toString(), gene, function);
        }
        return ImmutableKnownFusions.builder().gene("").eventType("").source("").sourceLink("").build();
    }

    @NotNull
    public static Signatures checkSignatures(@NotNull ViccEntry viccEntry, @NotNull String event) {
        if (SIGNATURE_MSI.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("MSI");
            return SignaturesExtractor.determineSignatures(viccEntry.source(), typeEvent.toString());
        } else if (SIGNATURE_HRD.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("HRD");
            return SignaturesExtractor.determineSignatures(viccEntry.source(), typeEvent.toString());
        } else if (SIGNATURE_MTL.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("MTL");
            return SignaturesExtractor.determineSignatures(viccEntry.source(), typeEvent.toString());
        } else if (SIGNATURE_MTB.contains(event)) {
            GenomicEvents typeEvent = GenomicEvents.genomicEvents("MTB");
            return SignaturesExtractor.determineSignatures(viccEntry.source(), typeEvent.toString());
        }

        return ImmutableSignatures.builder().eventType(Strings.EMPTY).source(Strings.EMPTY).sourceLink(Strings.EMPTY).build();
    }
}
