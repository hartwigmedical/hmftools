package com.hartwig.hmftools.knowledgebasegenerator.cnv;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.knowledgebasegenerator.eventtype.EventType;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneratingCNV {

    private static final Logger LOGGER = LogManager.getLogger(GeneratingCNV.class);

    private static final List<String> AMPLIFICATION =
            Lists.newArrayList("amplification", "Amplification", "Gain", "overexpression", "amp", "over exp");
    private static final List<String> DELETION =
            Lists.newArrayList("deletion", "Deletion", "Copy Number Loss", "Loss", "loss", "undexpression");

    public static void generatingCNVs(@NotNull ViccEntry viccEntries, @NotNull EventType type) {

        if (AMPLIFICATION.contains(type.eventType())) {
            // extract evidence items

        } else if (DELETION.contains(type.eventType())) {
            // extract evidence items
        }

    }
}

