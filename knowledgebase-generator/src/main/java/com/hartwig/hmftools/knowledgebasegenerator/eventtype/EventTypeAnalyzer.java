package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import java.util.List;
import java.util.regex.Pattern;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.brca.Brca;
import com.hartwig.hmftools.vicc.datamodel.cgi.Cgi;
import com.hartwig.hmftools.vicc.datamodel.civic.Civic;
import com.hartwig.hmftools.vicc.datamodel.jax.Jax;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKb;
import com.hartwig.hmftools.vicc.datamodel.pmkb.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.sage.Sage;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class EventTypeAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeAnalyzer.class);

    @NotNull
    public static List<EventType> determineEventType(@NotNull ViccEntry viccEntry) {
        KbSpecificObject kbSpecificObject = viccEntry.KbSpecificObject();
        String event = Strings.EMPTY;
        List<EventType> eventType = Lists.newArrayList();
        String sourceKnowledgebase = viccEntry.source();
        for (Feature feature : viccEntry.features()) {
            switch (sourceKnowledgebase) {
                case "oncokb": // extract info oncokb
                    OncoKb kbOncoKb = (OncoKb) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event.equals("NA")) {
                        event = feature.description().substring(feature.description().lastIndexOf(" ") + 1);
                        if (Pattern.compile("[0-9]").matcher(event).find()) {
                            LOGGER.info("manual curated: " + event + " to mutation");
                            event = "mutation";
                        }
                    }
                    LOGGER.info("event oncokb: " + event);
                    break;
                case "cgi": // extract info for cgi
                    Cgi kbCgi = (Cgi) kbSpecificObject;
                    event = feature.biomarkerType();
                    LOGGER.info("event cgi: " + event);
                    break;
                case "brca": // extract info for brca
                    Brca kbBrca = (Brca) kbSpecificObject;
                    event = Strings.EMPTY;
                    // TODO: extract event type
                    break;
                case "civic": // extract info for civic
                    Civic kbCivic = (Civic) kbSpecificObject;
                    event = feature.biomarkerType();
                    LOGGER.info("event civic: " + event);
                    // TODO: extract event type
                    break;
                case "jax": // extract info for jax
                    Jax kbJax = (Jax) kbSpecificObject;
                    event = Strings.EMPTY;
                    // TODO: extract event type
                    break;
                case "jax_trials": // extract info for jax trials
                    JaxTrials kbJaxTrials = (JaxTrials) kbSpecificObject;
                    event = Strings.EMPTY;
                    // TODO: extract event type
                    break;
                case "molecularmatch": // extract info for molecular match
                    MolecularMatch kbMolecularMatch = (MolecularMatch) kbSpecificObject;
                    event = Strings.EMPTY;
                    // TODO: extract event type
                    break;
                case "molecularmatch_trials": // extract info for molecular match trials
                    MolecularMatchTrials kbMolecularMatchTrials = (MolecularMatchTrials) kbSpecificObject;
                    event = Strings.EMPTY;
                    // TODO: extract event type
                    break;
                case "pmkb": // extract info for pmkb
                    Pmkb kbPmkb = (Pmkb) kbSpecificObject;
                    event = Strings.EMPTY;
                    // TODO: extract event type
                    break;
                case "sage": // extract info for sage
                    Sage kbSage = (Sage) kbSpecificObject;
                    event = Strings.EMPTY;
                    // TODO: extract event type
                    break;
                default:
                    LOGGER.info("Unknown knowledgebase");
            }
            if (event == null) {
                event = Strings.EMPTY;
            }
            eventType.add(ImmutableEventType.builder().gene("").eventType(event).build());
        }
        return eventType;
    }
}
