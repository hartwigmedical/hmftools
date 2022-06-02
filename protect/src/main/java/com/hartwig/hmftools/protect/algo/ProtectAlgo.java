package com.hartwig.hmftools.protect.algo;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordDataLoader;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDataLoader;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleDataLoader;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.protect.ProtectConfig;
import com.hartwig.hmftools.protect.evidence.ChordEvidence;
import com.hartwig.hmftools.protect.evidence.CopyNumberEvidence;
import com.hartwig.hmftools.protect.evidence.DisruptionEvidence;
import com.hartwig.hmftools.protect.evidence.FusionEvidence;
import com.hartwig.hmftools.protect.evidence.PersonalizedEvidenceFactory;
import com.hartwig.hmftools.protect.evidence.PurpleSignatureEvidence;
import com.hartwig.hmftools.protect.evidence.VariantEvidence;
import com.hartwig.hmftools.protect.evidence.VirusEvidence;
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
    private final VirusEvidence virusEvidenceFactory;
    @NotNull
    private final ChordEvidence chordEvidenceFactory;

    @NotNull
    public static ProtectAlgo build(@NotNull ActionableEvents actionableEvents, @NotNull Set<String> patientTumorDoids) {
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
                new PurpleSignatureEvidence(personalizedEvidenceFactory, actionableEvents.characteristics());
        VirusEvidence virusEvidenceFactory = new VirusEvidence(personalizedEvidenceFactory, actionableEvents.characteristics());
        ChordEvidence chordEvidenceFactory = new ChordEvidence(personalizedEvidenceFactory, actionableEvents.characteristics());

        return new ProtectAlgo(variantEvidenceFactory,
                copyNumberEvidenceFactory,
                disruptionEvidenceFactory,
                fusionEvidenceFactory,
                purpleSignatureEvidenceFactory,
                virusEvidenceFactory,
                chordEvidenceFactory);
    }

    private ProtectAlgo(@NotNull final VariantEvidence variantEvidenceFactory, @NotNull final CopyNumberEvidence copyNumberEvidenceFactory,
            @NotNull final DisruptionEvidence disruptionEvidenceFactory, @NotNull final FusionEvidence fusionEvidenceFactory,
            @NotNull final PurpleSignatureEvidence purpleSignatureEvidenceFactory, @NotNull final VirusEvidence virusEvidenceFactory,
            @NotNull final ChordEvidence chordEvidenceFactory) {
        this.variantEvidenceFactory = variantEvidenceFactory;
        this.copyNumberEvidenceFactory = copyNumberEvidenceFactory;
        this.disruptionEvidenceFactory = disruptionEvidenceFactory;
        this.fusionEvidenceFactory = fusionEvidenceFactory;
        this.purpleSignatureEvidenceFactory = purpleSignatureEvidenceFactory;
        this.virusEvidenceFactory = virusEvidenceFactory;
        this.chordEvidenceFactory = chordEvidenceFactory;
    }

    @NotNull
    public List<ProtectEvidence> run(@NotNull ProtectConfig config) throws IOException {
        PurpleData purpleData = loadPurpleData(config);
        LinxData linxData = loadLinxData(config);
        VirusInterpreterData virusInterpreterData = loadVirusInterpreterData(config);
        ChordAnalysis chordAnalysis = ChordDataLoader.load(config.chordPredictionTxt());

        return determineEvidence(purpleData, linxData, virusInterpreterData, chordAnalysis);
    }

    @NotNull
    private static PurpleData loadPurpleData(@NotNull ProtectConfig config) throws IOException {
        return PurpleDataLoader.load(config.tumorSampleId(),
                config.referenceSampleId(),
                null,
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleSomaticDriverCatalogTsv(),
                config.purpleSomaticVariantVcf(),
                config.purpleGermlineDriverCatalogTsv(),
                config.purpleGermlineVariantVcf(),
                config.purpleGeneCopyNumberTsv(),
                null,
                null,
                null);
    }

    @NotNull
    private static LinxData loadLinxData(@NotNull ProtectConfig config) throws IOException {
        return LinxDataLoader.load(config.linxFusionTsv(), config.linxBreakendTsv(), config.linxDriverCatalogTsv());
    }

    @NotNull
    private static VirusInterpreterData loadVirusInterpreterData(@NotNull ProtectConfig config) throws IOException {
        return VirusInterpreterDataLoader.load(config.annotatedVirusTsv());
    }

    @NotNull
    private List<ProtectEvidence> determineEvidence(@NotNull PurpleData purpleData, @NotNull LinxData linxData,
            @NotNull VirusInterpreterData virusInterpreterData, @NotNull ChordAnalysis chordAnalysis) {
        LOGGER.info("Evidence extraction started");
        List<ProtectEvidence> variantEvidence = variantEvidenceFactory.evidence(purpleData.reportableGermlineVariants(),
                purpleData.reportableSomaticVariants(),
                purpleData.unreportedSomaticVariants());
        printExtraction("somatic and germline variants", variantEvidence);
        List<ProtectEvidence> copyNumberEvidence =
                copyNumberEvidenceFactory.evidence(purpleData.reportableSomaticGainsLosses(), purpleData.unreportedSomaticGainsLosses());
        printExtraction("amplifications and deletions", copyNumberEvidence);
        List<ProtectEvidence> disruptionEvidence = disruptionEvidenceFactory.evidence(linxData.homozygousDisruptions());
        printExtraction("homozygous disruptions", disruptionEvidence);
        List<ProtectEvidence> fusionEvidence = fusionEvidenceFactory.evidence(linxData.reportableFusions(), linxData.unreportedFusions());
        printExtraction("fusions", fusionEvidence);
        List<ProtectEvidence> purpleSignatureEvidence = purpleSignatureEvidenceFactory.evidence(purpleData);
        printExtraction("purple signatures", purpleSignatureEvidence);
        List<ProtectEvidence> virusEvidence = virusEvidenceFactory.evidence(virusInterpreterData);
        printExtraction("viruses", virusEvidence);
        List<ProtectEvidence> chordEvidence = chordEvidenceFactory.evidence(chordAnalysis);
        printExtraction("chord", chordEvidence);

        List<ProtectEvidence> result = Lists.newArrayList();
        result.addAll(variantEvidence);
        result.addAll(copyNumberEvidence);
        result.addAll(disruptionEvidence);
        result.addAll(fusionEvidence);
        result.addAll(purpleSignatureEvidence);
        result.addAll(virusEvidence);
        result.addAll(chordEvidence);

        List<ProtectEvidence> consolidated = EvidenceConsolidation.consolidate(result);
        LOGGER.debug("Consolidated {} evidence items to {} unique evidence items", result.size(), consolidated.size());

        List<ProtectEvidence> reported = EvidenceReportingFunctions.applyReportingAlgo(consolidated);
        LOGGER.debug("Reduced reported evidence from {} items to {} items after applying reporting algo",
                reportedCount(consolidated),
                reportedCount(reported));

        List<ProtectEvidence> updatedForTrials = EvidenceReportingFunctions.reportOnLabelTrialsOnly(reported);
        LOGGER.debug("Reduced reported evidence from {} items to {} items by removing off-label trials",
                reportedCount(reported),
                reportedCount(updatedForTrials));

        List<ProtectEvidence> updatedForBlacklist = EvidenceReportingCuration.applyReportingBlacklist(updatedForTrials);
        LOGGER.debug("Reduced reported evidence from {} items to {} items by blacklisting specific evidence for reporting",
                reportedCount(updatedForTrials),
                reportedCount(updatedForBlacklist));

        return updatedForBlacklist;
    }

    private static void printExtraction(@NotNull String title, @NotNull List<ProtectEvidence> evidences) {
        Set<EvidenceKey> keys = EvidenceKey.buildKeySet(evidences);
        LOGGER.debug("Extracted {} evidence items for {} having {} keys ", evidences.size(), title, keys.size());
        for (EvidenceKey key : keys) {
            int count = evidences.stream().filter(x -> EvidenceKey.create(x).equals(key)).collect(Collectors.toList()).size();
            LOGGER.debug(" Resolved {} items for '{}'", count, key);
        }
    }

    private static int reportedCount(@NotNull List<ProtectEvidence> evidences) {
        return evidences.stream().filter(x -> x.reported()).collect(Collectors.toList()).size();
    }
}
