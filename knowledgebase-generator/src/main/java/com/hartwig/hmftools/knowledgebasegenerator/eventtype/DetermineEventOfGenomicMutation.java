package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.knowledgebasegenerator.actionability.gene.ActionableGene;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.CnvExtractor;
import com.hartwig.hmftools.knowledgebasegenerator.hotspot.HotspotExtractor;
import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class DetermineEventOfGenomicMutation {
    private static final Logger LOGGER = LogManager.getLogger(DetermineEventOfGenomicMutation.class);

    private static final List<String> AMPLIFICATION =
            Lists.newArrayList("Amplification", "Overexpression", "amp", "OVEREXPRESSION", "Transcript Amplification");
    private static final List<String> DELETION = Lists.newArrayList("Copy Number Loss", "Deletion", "del", "DELETION", "UNDEREXPRESSION");
    private static final Set<String> VARIANTS =
            Sets.newHashSet("missense_variant", "inframe_deletion", "inframe_insertion");
    private static final Set<String> FUSIONS = Sets.newHashSet();

    public static void checkGenomicEvent(@NotNull ViccEntry viccEntry, @NotNull EventType type, @NotNull HotspotExtractor hotspotExtractor)
            throws IOException, InterruptedException {
        //hotspotExtractor.extractHotspots(viccEntry);


        Source source = Source.sourceFromKnowledgebase(viccEntry.source());
        KbSpecificObject kbSpecificObject = viccEntry.KbSpecificObject();
        String gene = type.gene();
        String typeEvent = Strings.EMPTY;

        if (AMPLIFICATION.contains(type.eventType())) {
            typeEvent = "Amplification";
            ActionableGene actionableAmplification = CnvExtractor.determineInfoOfEvent(source, typeEvent, kbSpecificObject, gene, viccEntry);
            ActionableGene knownAmplification = CnvExtractor.determineInfoOfEvent(source, typeEvent, kbSpecificObject, gene, viccEntry);
        } else if (DELETION.contains(type.eventType())) {
            typeEvent = "Deletion";
            ActionableGene actionableDeletion = CnvExtractor.determineInfoOfEvent(source, typeEvent, kbSpecificObject, gene, viccEntry);
            ActionableGene knownDeletions = CnvExtractor.determineInfoOfEvent(source, typeEvent, kbSpecificObject, gene, viccEntry);
        } else if (VARIANTS.contains(type.eventType())) {
            // TODO: Determine hotspots
            hotspotExtractor.extractHotspots(viccEntry);
        } else if (FUSIONS.contains(type.eventType())) {
            // TODO: Determine fusions
        } else {
            LOGGER.info("skipping");
        }
    }
}
