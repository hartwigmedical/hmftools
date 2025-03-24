package com.hartwig.hmftools.common.doid;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.immutables.value.internal.$processor$.meta.$GsonMirrors;
import org.junit.Ignore;
import org.junit.Test;

public class DiseaseOntologyTest
{
    private static final String DOID_FILE_JSON = Resources.getResource("doid/example_doid.json").getPath();

    @Test
    public void canExtractDoidFromUrl()
    {
        String url = "http://purl.obolibrary.org/obo/DOID_345";
        assertEquals("345", DiseaseOntology.extractDoid(url));
    }

    @Ignore
    @Test
    public void canLoadDoidJsonFile() throws IOException
    {
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(DOID_FILE_JSON);

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

        assertTrue(doidEntry.meta().basicPropertyValues().isEmpty());
        assertTrue(doidEntry.logicalDefinitionAxioms().isEmpty());

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
        assertEquals("type", doidSynonym1.synonymType());

        DoidBasicPropertyValue doidBasicPropertyValue1 = doidNode1.doidMetadata().basicPropertyValues().get(0);
        assertEquals("http://www.geneontology.org/formats/oboInOwl#hasAlternativeId", doidBasicPropertyValue1.pred());
        assertEquals("DOID:8965", doidBasicPropertyValue1.val());

        DoidBasicPropertyValue doidBasicPropertyValue2 = doidNode1.doidMetadata().basicPropertyValues().get(1);
        assertEquals("http://www.geneontology.org/formats/oboInOwl#hasOBONamespace", doidBasicPropertyValue2.pred());
        assertEquals("disease_ontology", doidBasicPropertyValue2.val());

        assertTrue(doidNode1.doidMetadata().deprecated());
        assertEquals(Lists.newArrayList(), doidNode1.doidMetadata().comments());

        DoidNode doidNode2 = doidNodes.get(1);
        assertEquals("8717", doidNode2.doid());
        assertEquals("http://purl.obolibrary.org/obo/DOID_8717", doidNode2.url());
        assertEquals("decubitus ulcer", doidNode2.doidTerm());
        assertEquals("CLASS", doidNode2.type());

        DoidDefinition doidDefinition2 = doidNode2.doidMetadata().doidDefinition();
        assertEquals("Decubitus ulcer is a chronic ulcer of skin where the ulcer is an ulceration of "
                + "tissue deprived of adequate blood supply by prolonged pressure.", doidDefinition2.definitionVal());
        assertEquals(Lists.newArrayList(),
                doidDefinition2.definitionXrefs());

        List<String> subset2 = doidNode2.doidMetadata().subsets();
        assertEquals(Lists.newArrayList("http://purl.obolibrary.org/obo/doid#NCIthesaurus"), subset2);

        List<DoidXref> doidXrefs = doidNode2.doidMetadata().xrefs();
        assertEquals(ImmutableDoidXref.builder().val("NCI:C50706").build(), doidXrefs.get(0));
        assertEquals(ImmutableDoidXref.builder().val("MESH:D003668").build(), doidXrefs.get(1));
        assertEquals(ImmutableDoidXref.builder().val("ICD9CM:707.0").build(), doidXrefs.get(2));
        assertEquals(ImmutableDoidXref.builder().val("SNOMEDCT_US_2021_09_01:28103007").build(), doidXrefs.get(3));
        assertEquals(ImmutableDoidXref.builder().val("UMLS_CUI:C0011127").build(), doidXrefs.get(4));
        assertEquals(ImmutableDoidXref.builder().val("ICD10CM:L89").build(), doidXrefs.get(5));

        DoidSynonym doid2Synonym1 = doidNode2.doidMetadata().synonyms().get(0);
        assertEquals("hasExactSynonym", doid2Synonym1.pred());
        assertEquals("Decubitus ulcer any site", doid2Synonym1.val());
        assertNull(doid2Synonym1.synonymType());

        DoidSynonym doid2Synonym2 = doidNode2.doidMetadata().synonyms().get(1);
        assertEquals("hasExactSynonym", doid2Synonym2.pred());
        assertEquals("pressure ulcer", doid2Synonym2.val());

        DoidSynonym doid2Synonym3 = doidNode2.doidMetadata().synonyms().get(2);
        assertEquals("hasExactSynonym", doid2Synonym3.pred());
        assertEquals("pressure sores", doid2Synonym3.val());

        DoidSynonym doid2Synonym4 = doidNode2.doidMetadata().synonyms().get(3);
        assertEquals("hasExactSynonym", doid2Synonym4.pred());
        assertEquals("Decubitus (pressure) ulcer", doid2Synonym4.val());

        DoidSynonym doid2Synonym5 = doidNode2.doidMetadata().synonyms().get(4);
        assertEquals("hasRelatedSynonym", doid2Synonym5.pred());
        assertEquals("bedsore", doid2Synonym5.val());

        DoidBasicPropertyValue doid2BasicPropertyValue1 = doidNode2.doidMetadata().basicPropertyValues().get(0);
        assertEquals("http://www.geneontology.org/formats/oboInOwl#hasAlternativeId", doid2BasicPropertyValue1.pred());
        assertEquals("DOID:8808", doid2BasicPropertyValue1.val());

        DoidBasicPropertyValue doid2BasicPropertyValue2 = doidNode2.doidMetadata().basicPropertyValues().get(1);
        assertEquals("http://www.geneontology.org/formats/oboInOwl#hasAlternativeId", doid2BasicPropertyValue2.pred());
        assertEquals("DOID:9129", doid2BasicPropertyValue2.val());

        DoidBasicPropertyValue doid2BasicPropertyValue3 = doidNode2.doidMetadata().basicPropertyValues().get(2);
        assertEquals("http://www.geneontology.org/formats/oboInOwl#hasOBONamespace", doid2BasicPropertyValue3.pred());
        assertEquals("disease_ontology", doid2BasicPropertyValue3.val());

        DoidBasicPropertyValue doid2BasicPropertyValue4 = doidNode2.doidMetadata().basicPropertyValues().get(3);
        assertEquals("http://www.geneontology.org/formats/oboInOwl#hasAlternativeId", doid2BasicPropertyValue4.pred());
        assertEquals("DOID:9029", doid2BasicPropertyValue4.val());

        DoidBasicPropertyValue doid2BasicPropertyValue5 = doidNode2.doidMetadata().basicPropertyValues().get(4);
        assertEquals("http://www.geneontology.org/formats/oboInOwl#hasAlternativeId", doid2BasicPropertyValue5.pred());
        assertEquals("DOID:9002", doid2BasicPropertyValue5.val());

        assertNull(doidNode2.doidMetadata().deprecated());
        assertEquals(Lists.newArrayList("Xref MGI.\nOMIM mapping confirmed by DO. [SN]."), doidNode2.doidMetadata().comments());
    }
}