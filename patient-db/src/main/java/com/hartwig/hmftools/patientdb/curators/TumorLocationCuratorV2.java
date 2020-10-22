package com.hartwig.hmftools.patientdb.curators;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocationV2;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocationV2;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.semanticweb.HermiT.Reasoner;
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.io.FileDocumentSource;
import org.semanticweb.owlapi.model.*;

import java.util.logging.Level;


import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.semanticweb.owlapi.reasoner.InferenceType;
import org.semanticweb.owlapi.reasoner.OWLReasoner;

public class TumorLocationCuratorV2 implements CleanableCurator {

    private static final Logger LOGGER = LogManager.getLogger(BiopsySiteCurator.class);
    private static final IRI DISEASE_IRI = IRI.create("http://purl.obolibrary.org/obo/DOID_4");

    private static final String FIELD_DELIMITER = "\t";
    private static final String DOID_DELIMITER = ";";

    @NotNull
    private final Map<String, CuratedTumorLocationV2> tumorLocationMap = Maps.newHashMap();
    @NotNull
    private final Set<String> unusedSearchTerms;

    private OWLOntology createOntology(@NotNull String doidFile) throws OWLOntologyCreationException {
        java.util.logging.Logger logger = java.util.logging.Logger.getLogger("org.obolibrary.oboformat.parser.OBOFormatParser");
        logger.setLevel(Level.SEVERE);
        OWLOntologyManager ontologyManager = OWLManager.createOWLOntologyManager();
        OWLOntologyLoaderConfiguration config = new OWLOntologyLoaderConfiguration().setFollowRedirects(true)
                .setMissingImportHandlingStrategy(MissingImportHandlingStrategy.SILENT);

        return ontologyManager.loadOntologyFromOntologyDocument(new FileDocumentSource(new File(doidFile)), config);
    }

    private OWLReasoner createReasoner(@NotNull OWLOntology ontology) {
        OWLReasoner reasoner = new Reasoner.ReasonerFactory().createReasoner(ontology);
        reasoner.precomputeInferences(InferenceType.CLASS_HIERARCHY);
        return reasoner;
    }

    private void createDiseaseMapping(@NotNull OWLOntology ontology, @NotNull OWLReasoner reasoner) {

        OWLClass diseaseClass = ontology.getOWLOntologyManager().getOWLDataFactory().getOWLClass(DISEASE_IRI);
        LOGGER.info(diseaseClass);
//        val diseases = setOf(diseaseClass) + subClasses(diseaseClass, reasoner)
//        val synonyms = diseases.flatMap { disease -> getSynonyms(disease, ontology).map { Pair(it, disease) } }
//        return diseases.associateBy { getLabel(it, ontology) } + synonyms.toMap()
    }
    


    private List<String> extractDoidTerms(@NotNull String doidFile, @Nullable List<String> doids) throws OWLOntologyCreationException {
        //TODO read doid File

        OWLOntology ontology = createOntology(doidFile);
        OWLReasoner reasoner = createReasoner(ontology);
        createDiseaseMapping(ontology, reasoner);

        List<String> doidTerms = Lists.newArrayList();
        if (doids == null) {
            return null;
        } else {
            for (String doid : doids) {
                // TODO: implement logica for extract doidTerms
                doidTerms.add("doidTerms" + DOID_DELIMITER);
            }
            doidTerms.add(doidTerms.toString().substring(0, doidTerms.toString().length() - 1));
            return doidTerms;
        }
    }

    public TumorLocationCuratorV2(@NotNull String tumorLocationV2MappingTSV, @NotNull String doidFile)
            throws IOException, OWLOntologyCreationException {
        List<String> lines = Files.readAllLines(new File(tumorLocationV2MappingTSV).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_DELIMITER);
            String searchTerm = parts[0];
            String primaryTumorLocation = parts[1];
            String primaryTumorSubLocation = parts.length > 2 ? parts[2] : Strings.EMPTY;
            String primaryTumorType = parts.length > 3 ? parts[3] : Strings.EMPTY;
            String primaryTumorSubType = parts.length > 4 ? parts[4] : Strings.EMPTY;
            String primaryTumorExtraDetails = parts.length > 5 ? parts[5] : Strings.EMPTY;
            List<String> doids = parts.length > 6 ? Lists.newArrayList(parts[6].split(DOID_DELIMITER)) : Lists.newArrayList();
            tumorLocationMap.put(searchTerm,
                    ImmutableCuratedTumorLocationV2.builder()
                            .primaryTumorLocation(primaryTumorLocation)
                            .primaryTumorSubLocation(primaryTumorSubLocation)
                            .primaryTumorType(primaryTumorType)
                            .primaryTumorSubType(primaryTumorSubType)
                            .primaryTumorExtraDetails(primaryTumorExtraDetails)
                            .doids(doids)
                            .doidTerms(extractDoidTerms(doidFile, doids))
                            .searchTerm(searchTerm)
                            .build());
        }

        // Need to create a copy of the key set so that we can remove elements from it without affecting the curation.
        unusedSearchTerms = Sets.newHashSet(tumorLocationMap.keySet());
    }

    @NotNull
    public CuratedTumorLocationV2 search(@Nullable String searchTerm) {
        if (searchTerm != null && !searchTerm.isEmpty()) {
            unusedSearchTerms.remove(searchTerm);
            CuratedTumorLocationV2 result = tumorLocationMap.get(searchTerm);

            if (result != null) {
                return result;
            }
        }

        return ImmutableCuratedTumorLocationV2.builder().searchTerm(searchTerm).build();
    }

    @NotNull
    @Override
    public Set<String> unusedSearchTerms() {
        return unusedSearchTerms;
    }
}
