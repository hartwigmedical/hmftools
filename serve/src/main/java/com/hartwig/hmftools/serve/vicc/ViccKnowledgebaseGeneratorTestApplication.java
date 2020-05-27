package com.hartwig.hmftools.serve.vicc;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.vicc.copynumber.ActionableAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.eventtype.DetermineEventOfGenomicMutation;
import com.hartwig.hmftools.serve.vicc.eventtype.EventType;
import com.hartwig.hmftools.serve.vicc.eventtype.EventTypeAnalyzer;
import com.hartwig.hmftools.serve.vicc.fusion.KnownFusions;
import com.hartwig.hmftools.serve.vicc.signatures.Signatures;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;
import com.hartwig.hmftools.vicc.selection.ImmutableViccQuerySelection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ViccKnowledgebaseGeneratorTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ViccKnowledgebaseGeneratorTestApplication.class);

    public static void main(String[] args) throws IOException {
        String viccJsonPath = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";

        ViccSource source = ViccSource.ONCOKB;
        LOGGER.info("Reading VICC json from {} with source '{}'", viccJsonPath, source);
        List<ViccEntry> viccEntries = ViccJsonReader.readSelection(viccJsonPath,
                ImmutableViccQuerySelection.builder().addSourcesToFilterOn(source).build());
        LOGGER.info("Read {} entries", viccEntries.size());

        //Lists of known genomic events
        List<KnownAmplificationDeletion> listKnownAmplification = Lists.newArrayList();
        List<KnownAmplificationDeletion> listKnownDeletion = Lists.newArrayList();
        List<String> listKnownVariants = Lists.newArrayList();
        List<String> listKnownRange = Lists.newArrayList();
        List<KnownFusions> listKnownFusionPairs = Lists.newArrayList();
        List<KnownFusions> listKnownFusionPromiscuous = Lists.newArrayList();
        List<Signatures> listSignatures = Lists.newArrayList();

        //Lists of actionable genomic events
        List<ActionableAmplificationDeletion> listActionableDeletion = Lists.newArrayList();
        List<ActionableAmplificationDeletion> listActionableAmplification = Lists.newArrayList();

        for (ViccEntry viccEntry : viccEntries) {

            List<EventType> eventType = EventTypeAnalyzer.determineEventTypes(viccEntry);

            for (EventType type : eventType) {
                //  LOGGER.info("gene: " + type.gene() + " name: " + type.name() + " eventMap: " + type.eventMap() + " source: " + type.source());
                // Generating known events
                //TODO: map every genomic event to one object
                //TODO: if combined event use single event for determine known events

                for (Map.Entry<String, List<String>> entryDB : type.eventMap().entrySet()) {
                    for (String event : entryDB.getValue()) {

                        listKnownAmplification.add(DetermineEventOfGenomicMutation.checkKnownAmplification(viccEntry,
                                entryDB.getKey(),
                                event));
                        listKnownDeletion.add(DetermineEventOfGenomicMutation.checkKnownDeletion(viccEntry, entryDB.getKey(), event));
                        DetermineEventOfGenomicMutation.checkVariants(viccEntry, entryDB.getKey(), event);
                        DetermineEventOfGenomicMutation.checkRange(viccEntry, entryDB.getKey(), event);
                        listKnownFusionPairs.add(DetermineEventOfGenomicMutation.checkFusionsPairs(viccEntry, entryDB.getKey(), event));
                        listKnownFusionPromiscuous.add(DetermineEventOfGenomicMutation.checkFusionPromiscuous(viccEntry,
                                entryDB.getKey(),
                                event));
                        listSignatures.add(DetermineEventOfGenomicMutation.checkSignatures(viccEntry, event));
                    }
                }
            }
        }

        List<KnownAmplificationDeletion> listAmpsFilter = Lists.newArrayList();
        List<KnownAmplificationDeletion> listDelsFIlter = Lists.newArrayList();
        Set<String> uniqueAmps = Sets.newHashSet();
        Set<String> uniqueDels = Sets.newHashSet();
        for (KnownAmplificationDeletion amps : listKnownAmplification) {
            if (!amps.eventType().isEmpty()) {
                listAmpsFilter.add(amps);
                uniqueAmps.add(amps.gene());
            }
        }

        List<String> sortedUniqueAmps = Lists.newArrayList(uniqueAmps);
        Collections.sort(sortedUniqueAmps);

        for (KnownAmplificationDeletion dels : listKnownDeletion) {
            if (!dels.eventType().isEmpty()) {
                listDelsFIlter.add(dels);
                uniqueDels.add(dels.gene());
            }
        }
        List<String> sortedUniqueDels = Lists.newArrayList(uniqueDels);
        Collections.sort(sortedUniqueDels);

        List<Signatures> listSignaturesFilter = Lists.newArrayList();
        for (Signatures signatures : listSignatures) {
            if (!signatures.eventType().isEmpty()) {
                listSignaturesFilter.add(signatures);
            }
        }

        List<KnownFusions> listKnownFusionsPairsFilter = Lists.newArrayList();
        List<String> promiscusThree = Lists.newArrayList();
        List<String> promiscusFive = Lists.newArrayList();

        for (KnownFusions knownPairFusions : listKnownFusionPairs) {
            if (!knownPairFusions.eventType().isEmpty()) {
                listKnownFusionsPairsFilter.add(knownPairFusions);

                if (knownPairFusions.gene().equals("TRAC-NKX2-1")) {
                    promiscusFive.add("TRAC");
                    promiscusThree.add("NKX2-1");
                } else if (knownPairFusions.gene().contains("-")) {
                    promiscusFive.add(knownPairFusions.gene().split("-")[0]);
                    promiscusThree.add(knownPairFusions.gene().split("-")[1]);
                }
            }
        }

        Map<String, Integer> countsPromiscuousFive = Maps.newHashMap();
        for (String five : promiscusFive) {
            if (countsPromiscuousFive.containsKey(five)) {
                int count = countsPromiscuousFive.get(five) + 1;
                countsPromiscuousFive.put(five, count);
            } else {
                countsPromiscuousFive.put(five, 1);
            }
        }

        Set<String> promiscuousFiveGenes = Sets.newHashSet();
        for (KnownFusions five : listKnownFusionPromiscuous) {
            if (!five.eventType().isEmpty()) {
                if (countsPromiscuousFive.containsKey(five.gene())) {
                    if (countsPromiscuousFive.get(five.gene()) >= 3 || countsPromiscuousFive.get(five.gene()) >= 1) {
                        promiscuousFiveGenes.add(five.gene());
                    }
                }
            }
        }

        Map<String, Integer> countsPromiscuousThree = Maps.newHashMap();
        for (String three : promiscusThree) {
            if (countsPromiscuousThree.containsKey(three)) {
                int count = countsPromiscuousThree.get(three) + 1;
                countsPromiscuousThree.put(three, count);
            } else {
                countsPromiscuousThree.put(three, 1);
            }
        }

        Set<String> promiscuousThreeGenes = Sets.newHashSet();
        for (KnownFusions three : listKnownFusionPromiscuous) {
            if (!three.eventType().isEmpty()) {
                if (countsPromiscuousThree.containsKey(three.gene())) {
                    if (countsPromiscuousThree.get(three.gene()) >= 3 || countsPromiscuousThree.get(three.gene()) >= 1) {
                        promiscuousThreeGenes.add(three.gene());
                    }
                }
            }
        }

        List<ActionableAmplificationDeletion> listFilterActionableAmplifications = Lists.newArrayList();

        // If drug info/tumor location is known then variant is an actionable variant
        for (ActionableAmplificationDeletion actionableAmplification : listActionableAmplification) {
            if (actionableAmplification.level() != null && actionableAmplification.drug() != null
                    && actionableAmplification.drugType() != null && actionableAmplification.direction() != null
                    && actionableAmplification.sourceLink() != null && actionableAmplification.cancerType() != null) {
                listFilterActionableAmplifications.add(actionableAmplification);
            }
        }

        List<ActionableAmplificationDeletion> listFilterActionableDeletion = Lists.newArrayList();

        for (ActionableAmplificationDeletion actionableDeletion : listActionableDeletion) {
            if (actionableDeletion.level() != null && actionableDeletion.drug() != null && actionableDeletion.drugType() != null
                    && actionableDeletion.direction() != null && actionableDeletion.sourceLink() != null
                    && actionableDeletion.cancerType() != null) {
                listFilterActionableDeletion.add(actionableDeletion);
            }
        }

    }
}
