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
        String gene = Strings.EMPTY;
        String name = Strings.EMPTY;
        String description = Strings.EMPTY;
        Source source = Source.sourceFromKnowledgebase(viccEntry.source());
        Map<String, List<String>> eventMap = Maps.newHashMap();
        String eventInfo = Strings.EMPTY;

        List<EventType> eventType = Lists.newArrayList();

        for (Feature feature : viccEntry.features()) {
            switch (source) {
                case ONCOKB: // extract info oncokb
                    name = feature.name();
                    gene = feature.geneSymbol();
                    biomarkerType = feature.biomarkerType();
                    description = feature.description();

                    if (name.contains("Fusion")) {

                        if (name.split(" ").length == 2) {
                            gene = name.split(" ")[0];
                            name = name.split(" ")[1];
                        }
                    }

                    eventMap.put(gene, Lists.newArrayList(name));

                    if (name.isEmpty()){
                        LOGGER.warn(
                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' on source '{}' ",
                                feature.name(),
                                gene,
                                biomarkerType, source);
                    }

                    break;
                case CGI: // extract info for cgi
                    name = feature.name();
                    biomarkerType = feature.biomarkerType();
                    description = feature.description();

                    if (name.contains("+")) {
                        String[] combinedEventConvertToSingleEvent = name.split(" \\+ ", 2);
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
                        if (name.contains(":")) {
                            gene = name.split(":", 2)[0];
                            eventInfo = name.split(":", 2)[1];
                        } else {
                            gene = name.split(" ", 2)[0];
                            if (name.split(" ", 2).length == 1) {
                                eventInfo = "Fusion";
                            } else {
                                eventInfo = name.split(" ", 2)[1];
                            }
                        }
                        eventMap.put(gene, Lists.newArrayList(eventInfo));
                    }
                    if (eventInfo.isEmpty()){
                        LOGGER.warn(
                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' on source '{}'",
                                feature.name(),
                                gene,
                                biomarkerType, source);
                    }

                    break;
                case CIVIC: // extract info for civic
                    name = feature.name();
                    biomarkerType = feature.biomarkerType() == null ? "null" : feature.biomarkerType();
                    gene = feature.geneSymbol();
                    description = feature.description();

                    eventMap.put(gene, Lists.newArrayList(name));

                    if (name.isEmpty()){
                        LOGGER.warn(
                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' on source '{}'",
                                feature.name(),
                                gene,
                                biomarkerType, source);
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
                    .combinedEvent(combinedEvent)
                    .biomarkerType(biomarkerType)
                    .gene(gene)
                    .name(name)
                    .description(description)
                    .source(source)
                    .eventMap(eventMap)
                    .event(Strings.EMPTY)
                    .interpretEventType(Lists.newArrayList()).build());
        }
        return eventType;
    }
}
