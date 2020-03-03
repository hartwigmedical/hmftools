package com.hartwig.hmftools.knowledgebasegenerator.cnv;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.knowledgebasegenerator.eventtype.EventType;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CnvExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CnvExtractor.class);

    private static final List<String> AMPLIFICATION =
            Lists.newArrayList("Amplification", "Overexpression", "amp", "OVEREXPRESSION", "Transcript Amplification");
    private static final List<String> DELETION =
            Lists.newArrayList("Copy Number Loss", "Deletion", "del", "DELETION", "UNDEREXPRESSION");

    public static void extractingCNVs(@NotNull ViccEntry viccEntries, @NotNull EventType type) {

        if (AMPLIFICATION.contains(type.eventType())) {
            // extract evidence items

        } else if (DELETION.contains(type.eventType())) {
            // extract evidence items
        }

    }
}

