package com.hartwig.hmftools.protect;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.protect.bachelor.BachelorData;
import com.hartwig.hmftools.protect.evidence.ChordEvidence;
import com.hartwig.hmftools.protect.evidence.CopyNumberEvidence;
import com.hartwig.hmftools.protect.evidence.DisruptionEvidence;
import com.hartwig.hmftools.protect.evidence.FusionEvidence;
import com.hartwig.hmftools.protect.evidence.PersonalizedEvidenceFactory;
import com.hartwig.hmftools.protect.evidence.ProtectEvidenceFunctions;
import com.hartwig.hmftools.protect.evidence.PurpleSignatureEvidence;
import com.hartwig.hmftools.protect.evidence.VariantEvidence;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.serve.actionability.ActionableEvents;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ProtectAlgo {

    private static final Logger LOGGER = LogManager.getLogger(ProtectAlgo.class);

    @NotNull
    private final VariantEvidence variantEvidenceFactory;
    @NotNull
    private final CopyNumberEvidence copyNumberEvidenceFactory;
    @NotNull
    private final DisruptionEvidence disruptionEvidenceFactory;
    @NotNull
    private final FusionEvidence fusionEvidenceFactory;
    @NotNull
    private final PurpleSignatureEvidence purpleSignatureEvidenceFactory;
    @NotNull
    private final ChordEvidence chordEvidenceFactory;

    @NotNull
    public static ProtectAlgo buildAlgoFromServeActionability(@NotNull ActionableEvents actionableEvents,
            @NotNull Set<String> patientTumorDoids) {
        PersonalizedEvidenceFactory personalizedEvidenceFactory = new PersonalizedEvidenceFactory(patientTumorDoids);

        VariantEvidence variantEvidenceFactory = new VariantEvidence(personalizedEvidenceFactory,
                actionableEvents.hotspots(),
                actionableEvents.ranges(),
                actionableEvents.genes());
        CopyNumberEvidence copyNumberEvidenceFactory = new CopyNumberEvidence(personalizedEvidenceFactory, actionableEvents.genes());
        DisruptionEvidence disruptionEvidenceFactory = new DisruptionEvidence(personalizedEvidenceFactory, actionableEvents.genes());
        FusionEvidence fusionEvidenceFactory =
                new FusionEvidence(personalizedEvidenceFactory, actionableEvents.genes(), actionableEvents.fusions());
        PurpleSignatureEvidence purpleSignatureEvidenceFactory =
                new PurpleSignatureEvidence(personalizedEvidenceFactory, actionableEvents.signatures());
        ChordEvidence chordEvidenceFactory = new ChordEvidence(personalizedEvidenceFactory, actionableEvents.signatures());

        return new ProtectAlgo(variantEvidenceFactory,
                copyNumberEvidenceFactory,
                disruptionEvidenceFactory,
                fusionEvidenceFactory,
                purpleSignatureEvidenceFactory,
                chordEvidenceFactory);
    }

    private ProtectAlgo(@NotNull final VariantEvidence variantEvidenceFactory, @NotNull final CopyNumberEvidence copyNumberEvidenceFactory,
            @NotNull final DisruptionEvidence disruptionEvidenceFactory, @NotNull final FusionEvidence fusionEvidenceFactory,
            @NotNull final PurpleSignatureEvidence purpleSignatureEvidenceFactory, @NotNull final ChordEvidence chordEvidenceFactory) {
        this.variantEvidenceFactory = variantEvidenceFactory;
        this.copyNumberEvidenceFactory = copyNumberEvidenceFactory;
        this.disruptionEvidenceFactory = disruptionEvidenceFactory;
        this.fusionEvidenceFactory = fusionEvidenceFactory;
        this.purpleSignatureEvidenceFactory = purpleSignatureEvidenceFactory;
        this.chordEvidenceFactory = chordEvidenceFactory;
    }

    @NotNull
    public List<ProtectEvidence> determineEvidence(@NotNull PurpleData purpleData, @NotNull LinxData linxData,
            @NotNull BachelorData bachelorData, @NotNull ChordAnalysis chordAnalysis) {
        List<ProtectEvidence> variantEvidence =
                variantEvidenceFactory.evidence(bachelorData.germlineVariants(), purpleData.somaticVariants());
        printExtraction("somatic and germline variants", variantEvidence);
        List<ProtectEvidence> copyNumberEvidence = copyNumberEvidenceFactory.evidence(purpleData.copyNumberAlterations());
        printExtraction("amplifications and deletions", copyNumberEvidence);
        List<ProtectEvidence> disruptionEvidence = disruptionEvidenceFactory.evidence(linxData.homozygousDisruptions());
        printExtraction("homozygous disruptions", disruptionEvidence);
        List<ProtectEvidence> fusionEvidence = fusionEvidenceFactory.evidence(linxData.fusions());
        printExtraction("fusions", fusionEvidence);
        List<ProtectEvidence> purpleSignatureEvidence = purpleSignatureEvidenceFactory.evidence(purpleData);
        printExtraction("purple signatures", purpleSignatureEvidence);
        List<ProtectEvidence> chordEvidence = chordEvidenceFactory.evidence(chordAnalysis);
        printExtraction("chord", chordEvidence);

        List<ProtectEvidence> result = Lists.newArrayList();
        result.addAll(variantEvidence);
        result.addAll(copyNumberEvidence);
        result.addAll(disruptionEvidence);
        result.addAll(fusionEvidence);
        result.addAll(purpleSignatureEvidence);
        result.addAll(chordEvidence);

        List<ProtectEvidence> consolidated = ProtectEvidenceFunctions.consolidate(result);
        LOGGER.debug(" Consolidated {} evidences to {} unique evidences", result.size(), consolidated.size());

        return ProtectEvidenceFunctions.reportHighest(consolidated);
    }

    private static void printExtraction(@NotNull String title, @NotNull List<ProtectEvidence> evidences) {
        Set<String> events = evidences.stream().map(x -> x.genomicEvent()).collect(Collectors.toSet());
        LOGGER.debug("Extracted {} evidences for {} based off {} genomic events", evidences.size(), title, events.size());
        for (String event : events) {
            int count = evidences.stream().filter(x -> x.genomicEvent().equals(event)).collect(Collectors.toList()).size();
            LOGGER.debug(" Resolved {} evidences for '{}'", count, event);
        }
    }
}
