package com.hartwig.hmftools.common.doid;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

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

        assertEquals(DiseaseOntology.ID_TO_READ, doidEntry.id());

        List<DoidEdge> doidEdges = doidEntry.edges();
        assertEquals(9, doidEdges.size());

        DoidEdge doidEdge1 = doidEdges.get(0);
        assertEquals("http://purl.obolibrary.org/obo/DOID_8717", doidEdge1.subject());
        assertEquals("is_a", doidEdge1.predicate());
        assertEquals("http://purl.obolibrary.org/obo/DOID_8549", doidEdge1.object());

        DoidEdge doidEdge2 = doidEdges.get(1);
        assertEquals("http://purl.obolibrary.org/obo/CHEBI_50906", doidEdge2.subject());
        assertEquals("is_a", doidEdge2.predicate());
        assertEquals("http://purl.obolibrary.org/obo/doid#chebi", doidEdge2.object());

        assertTrue(doidEntry.meta().subsets().isEmpty());
        assertTrue(doidEntry.meta().xrefs().isEmpty());
        assertTrue(doidEntry.meta().basicPropertyValues().isEmpty());
        assertTrue(doidEntry.equivalentNodesSets().isEmpty());
        assertTrue(doidEntry.logicalDefinitionAxioms().isEmpty());
        assertTrue(doidEntry.domainRangeAxioms().isEmpty());
        assertTrue(doidEntry.propertyChainAxioms().isEmpty());

        List<DoidNode> doidNodes = doidEntry.nodes();
        assertEquals(2, doidNodes.size());

        DoidNode doidNode1 = doidNodes.get(0);
        assertEquals("8718", doidNode1.doid());
        assertEquals("http://purl.obolibrary.org/obo/DOID_8718", doidNode1.url());
        assertEquals("obsolete carcinoma in situ of respiratory system", doidNode1.doidTerm());
        assertEquals("CLASS", doidNode1.type());

        DoidDefinition doidDefinition1 = doidNode1.doidMetadata().doidDefinition();
        assertEquals("A carcinoma in situ that is characterized by the spread of cancer in the respiratory "
                + "system and the lack of invasion of surrounding tissues.", doidDefinition1.definitionVal());
        assertEquals(Lists.newArrayList("url:http://en.wikipedia.org/wiki/Carcinoma_in_situ"), doidDefinition1.definitionXrefs());

        DoidSynonym doidSynonym1 = doidNode1.doidMetadata().synonyms().get(0);
        assertEquals("hasExactSynonym", doidSynonym1.pred());
        assertEquals("carcinoma in situ of respiratory tract (disorder)", doidSynonym1.val());
        assertTrue(doidSynonym1.xrefs().isEmpty());

        DoidBasicPropertyValue doidBasicPropertyValue1 = doidNode1.doidMetadata().basicPropertyValues().get(0);
        assertEquals("http://www.geneontology.org/formats/oboInOwl#hasAlternativeId", doidBasicPropertyValue1.pred());
        assertEquals("DOID:8965", doidBasicPropertyValue1.val());

        DoidBasicPropertyValue doidBasicPropertyValue2 = doidNode1.doidMetadata().basicPropertyValues().get(1);
        assertEquals("http://www.w3.org/2002/07/owl#deprecated", doidBasicPropertyValue2.pred());
        assertEquals("true", doidBasicPropertyValue2.val());

        DoidBasicPropertyValue doidBasicPropertyValue3 = doidNode1.doidMetadata().basicPropertyValues().get(2);
        assertEquals("http://www.geneontology.org/formats/oboInOwl#hasOBONamespace", doidBasicPropertyValue3.pred());
        assertEquals("disease_ontology", doidBasicPropertyValue3.val());

        DoidNode doidNode2 = doidNodes.get(1);
        assertEquals("8717", doidNode2.doid());
        assertEquals("http://purl.obolibrary.org/obo/DOID_8717", doidNode2.url());
        assertEquals("decubitus ulcer", doidNode2.doidTerm());
        assertEquals("CLASS", doidNode2.type());

        DoidDefinition doidDefinition2 = doidNode2.doidMetadata().doidDefinition();
        assertEquals("Decubitus ulcer is a chronic ulcer of skin where the ulcer is an ulceration of "
                + "tissue deprived of adequate blood supply by prolonged pressure.", doidDefinition2.definitionVal());
        assertEquals(Lists.newArrayList("url:http://www2.merriam-webster.com/cgi-bin/mwmednlm?book=Medical&va=bedsore"),
                doidDefinition2.definitionXrefs());

        List<String> subset2 = doidNode2.doidMetadata().subsets();
        assertEquals(subset2, Lists.newArrayList("http://purl.obolibrary.org/obo/doid#NCIthesaurus"));

        List<DoidXref> doidXrefs = doidNode2.doidMetadata().xrefs();
        assertEquals(doidXrefs.get(0), ImmutableDoidXref.builder().val("NCI:C50706").build());
        assertEquals(doidXrefs.get(1), ImmutableDoidXref.builder().val("MESH:D003668").build());
        assertEquals(doidXrefs.get(2), ImmutableDoidXref.builder().val("ICD9CM:707.0").build());
        assertEquals(doidXrefs.get(3), ImmutableDoidXref.builder().val("UMLS_CUI:C0011127").build());
        assertEquals(doidXrefs.get(4), ImmutableDoidXref.builder().val("SNOMEDCT_US_2020_03_01:28103007").build());
        assertEquals(doidXrefs.get(5), ImmutableDoidXref.builder().val("ICD10CM:L89").build());

        DoidSynonym doid2Synonym1 = doidNode2.doidMetadata().synonyms().get(0);
        assertEquals(doid2Synonym1.pred(), "hasExactSynonym");
        assertEquals(doid2Synonym1.val(), "Decubitus ulcer any site");
        assertEquals(doid2Synonym1.xrefs(), Lists.newArrayList());

        DoidSynonym doid2Synonym2 = doidNode2.doidMetadata().synonyms().get(1);
        assertEquals(doid2Synonym2.pred(), "hasExactSynonym");
        assertEquals(doid2Synonym2.val(), "pressure ulcer");
        assertEquals(doid2Synonym2.xrefs(), Lists.newArrayList());

        DoidSynonym doid2Synonym3 = doidNode2.doidMetadata().synonyms().get(2);
        assertEquals(doid2Synonym3.pred(), "hasExactSynonym");
        assertEquals(doid2Synonym3.val(), "pressure sores");
        assertEquals(doid2Synonym3.xrefs(), Lists.newArrayList());

        DoidSynonym doid2Synonym4 = doidNode2.doidMetadata().synonyms().get(3);
        assertEquals(doid2Synonym4.pred(), "hasExactSynonym");
        assertEquals(doid2Synonym4.val(), "Decubitus (pressure) ulcer");
        assertEquals(doid2Synonym4.xrefs(), Lists.newArrayList());

        DoidSynonym doid2Synonym5 = doidNode2.doidMetadata().synonyms().get(4);
        assertEquals(doid2Synonym5.pred(), "hasRelatedSynonym");
        assertEquals(doid2Synonym5.val(), "bedsore");
        assertEquals(doid2Synonym5.xrefs(), Lists.newArrayList());

        DoidBasicPropertyValue doid2BasicPropertyValue1 = doidNode2.doidMetadata().basicPropertyValues().get(0);
        assertEquals(doid2BasicPropertyValue1.pred(), "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId");
        assertEquals(doid2BasicPropertyValue1.val(), "DOID:8808");

        DoidBasicPropertyValue doid2BasicPropertyValue2 = doidNode2.doidMetadata().basicPropertyValues().get(1);
        assertEquals(doid2BasicPropertyValue2.pred(), "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId");
        assertEquals(doid2BasicPropertyValue2.val(), "DOID:9129");

        DoidBasicPropertyValue doid2BasicPropertyValue3 = doidNode2.doidMetadata().basicPropertyValues().get(2);
        assertEquals(doid2BasicPropertyValue3.pred(), "http://www.geneontology.org/formats/oboInOwl#hasOBONamespace");
        assertEquals(doid2BasicPropertyValue3.val(), "disease_ontology");

        DoidBasicPropertyValue doid2BasicPropertyValue4 = doidNode2.doidMetadata().basicPropertyValues().get(3);
        assertEquals(doid2BasicPropertyValue4.pred(), "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId");
        assertEquals(doid2BasicPropertyValue4.val(), "DOID:9029");

        DoidBasicPropertyValue doid2BasicPropertyValue5 = doidNode2.doidMetadata().basicPropertyValues().get(4);
        assertEquals(doid2BasicPropertyValue5.pred(), "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId");
        assertEquals(doid2BasicPropertyValue5.val(), "DOID:9002");
    }
}