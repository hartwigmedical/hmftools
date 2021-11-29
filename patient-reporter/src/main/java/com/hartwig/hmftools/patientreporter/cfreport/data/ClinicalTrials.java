package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ClinicalTrials {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalTrials.class);

    private ClinicalTrials() {
    }

    public static int uniqueEventCount(@NotNull List<ProtectEvidence> trials) {
        Set<String> events = Sets.newHashSet();
        for (ProtectEvidence trial : trials) {
            events.add(trial.genomicEvent());
        }
        return events.size();
    }

    public static int uniqueTrialCount(@NotNull List<ProtectEvidence> trials) {
        Set<String> acronyms = Sets.newHashSet();
        for (ProtectEvidence trial : trials) {
            acronyms.add(trial.treatment());
        }
        return acronyms.size();
    }

    @NotNull
    public static String createLinkiClusion(@NotNull ProtectEvidence evidence) {
        String link = Strings.EMPTY;
        for (String url : evidence.urls()) {
            if (url.contains("trial-eye")) {
                link = url;
            }
        }
        //We assume iClusion has one link
        return link;
    }
}