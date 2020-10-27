package com.hartwig.hmftools.common.doid;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class DiseaseOntologyTest {

    private static final String DOID_FILE_JSON = Resources.getResource("doid/example_doid.json").getPath();

    @Test
    public void canExtractDoidFromUrl() {
        String url = "http://purl.obolibrary.org/obo/DOID_345";
        assertEquals("345", DiseaseOntology.extractDoid(url));
    }

    @Test
    public void canLoadDoidJsonFile() throws IOException {

        DoidEntry doidEntry = DiseaseOntology.readDoidJsonFile(DOID_FILE_JSON);
        DoidGraphMetaData doidGraphMetaData = ImmutableDoidGraphMetaData.builder()
                .subsets(Lists.newArrayList())
                .xrefs(Lists.newArrayList())
                .basicPropertyValues(Lists.newArrayList())
                .build();
        List<DoidNode> doidNodes = doidEntry.doidNodes();
        DoidEdges doidEdges = doidEntry.edges();

        assertEquals("http://purl.obolibrary.org/obo/doid/obo/ext.owl", doidEntry.id());
        assertEquals(doidGraphMetaData, doidEntry.meta());
        assertEquals(Lists.newArrayList(), doidEntry.equivalentNodesSets());
        assertEquals(Lists.newArrayList(), doidEntry.logicalDefinitionAxioms());
        assertEquals(Lists.newArrayList(), doidEntry.domainRangeAxioms());
        assertEquals(Lists.newArrayList(), doidEntry.propertyChainAxioms());

        assertEquals(2, doidNodes.size());

        DoidNode doidEntry1 = doidNodes.get(0);
        assertEquals(doidEntry1.doid(), "8718");
        assertEquals(doidEntry1.url(), "http://purl.obolibrary.org/obo/DOID_8718");
        assertEquals(doidEntry1.doidTerm(), "obsolete carcinoma in situ of respiratory system");
        assertEquals(doidEntry1.type(), "CLASS");

        DoidDefinition doidDefinition1 = doidEntry1.doidMetadata().doidDefinition();
        assertEquals(doidDefinition1.definitionVal(),
                "A carcinoma in situ that is characterized by the spread of cancer in the respiratory "
                        + "system and the lack of invasion of surrounding tissues.");
        assertEquals(doidDefinition1.definitionXrefs(), Lists.newArrayList("url:http://en.wikipedia.org/wiki/Carcinoma_in_situ"));

        DoidSynonym doidSynonym1 = doidEntry1.doidMetadata().synonyms().get(0);
        assertEquals(doidSynonym1.pred(), "hasExactSynonym");
        assertEquals(doidSynonym1.val(), "carcinoma in situ of respiratory tract (disorder)");
        assertEquals(doidSynonym1.xrefs(), Lists.newArrayList());

        DoidBasicPropertyValue doidBasicPropertyValue1 = doidEntry1.doidMetadata().basicPropertyValues().get(0);
        assertEquals(doidBasicPropertyValue1.pred(), "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId");
        assertEquals(doidBasicPropertyValue1.val(), "DOID:8965");

        DoidBasicPropertyValue doidBasicPropertyValue2 = doidEntry1.doidMetadata().basicPropertyValues().get(1);
        assertEquals(doidBasicPropertyValue2.pred(), "http://www.w3.org/2002/07/owl#deprecated");
        assertEquals(doidBasicPropertyValue2.val(), "true");

        DoidBasicPropertyValue doidBasicPropertyValue3 = doidEntry1.doidMetadata().basicPropertyValues().get(2);
        assertEquals(doidBasicPropertyValue3.pred(), "http://www.geneontology.org/formats/oboInOwl#hasOBONamespace");
        assertEquals(doidBasicPropertyValue3.val(), "disease_ontology");

        DoidNode doidEntry2 = doidNodes.get(1);
        assertEquals(doidEntry2.url(), "http://purl.obolibrary.org/obo/DOID_8717");
        assertEquals(doidEntry2.doid(), "8717");
        assertEquals(doidEntry2.doidTerm(), "decubitus ulcer");
    }
}