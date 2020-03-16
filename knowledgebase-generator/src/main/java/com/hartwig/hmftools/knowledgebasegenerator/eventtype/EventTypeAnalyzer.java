package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class EventTypeAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeAnalyzer.class);

    @NotNull
    public static List<EventType> determineEventType(@NotNull ViccEntry viccEntry) {

        boolean combinedEvent = false;
        String biomarkerType = Strings.EMPTY;
        String eventInfo = Strings.EMPTY;
        String gene = Strings.EMPTY;
        String eventSource = Strings.EMPTY;
        List<String> event = Lists.newArrayList();
        Map<String, List<String>> eventMap = Maps.newHashMap();

        List<EventType> eventType = Lists.newArrayList();

        Source type = Source.sourceFromKnowledgebase(viccEntry.source());

        for (Feature feature : viccEntry.features()) {
            switch (type) {
                case ONCOKB: // extract info oncokb
                    eventSource = feature.name();
                    gene = feature.geneSymbol();
                    biomarkerType = feature.biomarkerType();

                    if (eventSource.contains("Fusion")) {

                        if (eventSource.split(" ").length == 2) {
                            gene = eventSource.split(" ")[0];
                            eventSource = eventSource.split(" ")[1];
                        }
                    }

                    eventMap.put(gene, Lists.newArrayList(eventSource));

                    if (eventSource.isEmpty()){
                        LOGGER.warn(
                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' on source '{}' ",
                                feature.name(),
                                gene,
                                biomarkerType, type);
                    }

                    break;
                case CGI: // extract info for cgi
                    eventSource = feature.name();
                    biomarkerType = feature.biomarkerType();
                    if (eventSource.contains("+")) {
                        String[] combinedEventConvertToSingleEvent = eventSource.split(" \\+ ", 2);
                        gene = combinedEventConvertToSingleEvent[0].split(" ", 2)[0];
                        eventInfo = combinedEventConvertToSingleEvent[0].split(" ", 2)[1];

                        String geneCombined = combinedEventConvertToSingleEvent[1].split(" ", 2)[0];
                        String eventInfoCombined = combinedEventConvertToSingleEvent[1].split(" ", 2)[1];

                        if (combinedEventConvertToSingleEvent.length == 2) {
                            combinedEvent = true;

                            if (eventMap.size() == 0) {
                                eventMap.put(gene, Lists.newArrayList(eventInfo));
                                if (eventMap.containsKey(geneCombined)) {
                                    eventMap.put(geneCombined,
                                            Lists.newArrayList(eventInfo, eventInfoCombined));
                                } else {
                                    eventMap.put(gene, Lists.newArrayList(eventInfo));
                                    eventMap.put(geneCombined, Lists.newArrayList(eventInfoCombined));
                                }
                            }
                        }
                    } else {
                        if (eventSource.contains(":")) {
                            gene = eventSource.split(":", 2)[0];
                            eventInfo = eventSource.split(":", 2)[1];
                        } else {
                            gene = eventSource.split(" ", 2)[0];
                            if (eventSource.split(" ", 2).length == 1) {
                                eventInfo = "Fusion";
                            } else {
                                eventInfo = eventSource.split(" ", 2)[1];
                            }
                        }
                        eventMap.put(gene, Lists.newArrayList(eventInfo));
                    }
                    if (eventInfo.isEmpty()){
                        LOGGER.warn(
                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' on source '{}'",
                                feature.name(),
                                gene,
                                biomarkerType, type);
                    }

                    break;
                case CIVIC: // extract info for civic
                    eventSource = feature.name();
                    biomarkerType = feature.biomarkerType() == null ? "null" : feature.biomarkerType();
                    gene = feature.geneSymbol();
                    eventMap.put(gene, Lists.newArrayList(eventSource));

                    if (eventSource.isEmpty()){
                        LOGGER.warn(
                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' on source '{}'",
                                feature.name(),
                                gene,
                                biomarkerType, type);
                    }

                    break;
                case BRCA: // extract info for brca
                    break;
                case JAX: // extract info for jax
                    break;
                case JAX_TRIALS: // extract info for jax trials
                    break;
                case MOLECULARMATCH: // extract info for molecular match
                    break;
                case MOLECULARMATCH_TRIALS: // extract info for molecular match trials
                    break;
                case PMKB: // extract info for pmkb
                    break;
                case SAGE: // extract info for sage
                    break;
                default:
                    LOGGER.warn("Unknown knowledgebase!");
            }

            eventType.add(ImmutableEventType.builder()
                    .gene(gene)
                    .combinedEvent(combinedEvent)
                    .event(eventSource)
                    .interpretEventType(Lists.newArrayList())
                    .biomarkerType(biomarkerType)
                    .build());
        }
        return eventType;
    }
}
