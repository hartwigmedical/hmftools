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
        assertEquals(doidEntry2.doid(), "8717");
        assertEquals(doidEntry2.url(), "http://purl.obolibrary.org/obo/DOID_8717");
        assertEquals(doidEntry2.doidTerm(), "decubitus ulcer");
        assertEquals(doidEntry2.type(), "CLASS");

        DoidDefinition doidDefinition2 = doidEntry2.doidMetadata().doidDefinition();
        assertEquals(doidDefinition2.definitionVal(),
                "Decubitus ulcer is a chronic ulcer of skin where the ulcer is an ulceration of "
                        + "tissue deprived of adequate blood supply by prolonged pressure.");
        assertEquals(doidDefinition2.definitionXrefs(),
                Lists.newArrayList("url:http://www2.merriam-webster.com/cgi-bin/mwmednlm?book=Medical&va=bedsore"));

        List<String> subset2 = doidEntry2.doidMetadata().subsets();
        assertEquals(subset2, Lists.newArrayList("http://purl.obolibrary.org/obo/doid#NCIthesaurus"));

        List<DoidXref> doidXrefs = doidEntry2.doidMetadata().xrefs();
        assertEquals(doidXrefs.get(0), ImmutableDoidXref.builder().val("NCI:C50706").build());
        assertEquals(doidXrefs.get(1), ImmutableDoidXref.builder().val("MESH:D003668").build());
        assertEquals(doidXrefs.get(2), ImmutableDoidXref.builder().val("ICD9CM:707.0").build());
        assertEquals(doidXrefs.get(3), ImmutableDoidXref.builder().val("UMLS_CUI:C0011127").build());
        assertEquals(doidXrefs.get(4), ImmutableDoidXref.builder().val("SNOMEDCT_US_2020_03_01:28103007").build());
        assertEquals(doidXrefs.get(5), ImmutableDoidXref.builder().val("ICD10CM:L89").build());

        DoidSynonym doid2Synonym1 = doidEntry2.doidMetadata().synonyms().get(0);
        assertEquals(doid2Synonym1.pred(), "hasExactSynonym");
        assertEquals(doid2Synonym1.val(), "Decubitus ulcer any site");
        assertEquals(doid2Synonym1.xrefs(), Lists.newArrayList());

        DoidSynonym doid2Synonym2 = doidEntry2.doidMetadata().synonyms().get(1);
        assertEquals(doid2Synonym2.pred(), "hasExactSynonym");
        assertEquals(doid2Synonym2.val(), "pressure ulcer");
        assertEquals(doid2Synonym2.xrefs(), Lists.newArrayList());

        DoidSynonym doid2Synonym3 = doidEntry2.doidMetadata().synonyms().get(2);
        assertEquals(doid2Synonym3.pred(), "hasExactSynonym");
        assertEquals(doid2Synonym3.val(), "pressure sores");
        assertEquals(doid2Synonym3.xrefs(), Lists.newArrayList());

        DoidSynonym doid2Synonym4 = doidEntry2.doidMetadata().synonyms().get(3);
        assertEquals(doid2Synonym4.pred(), "hasExactSynonym");
        assertEquals(doid2Synonym4.val(), "Decubitus (pressure) ulcer");
        assertEquals(doid2Synonym4.xrefs(), Lists.newArrayList());

        DoidSynonym doid2Synonym5 = doidEntry2.doidMetadata().synonyms().get(4);
        assertEquals(doid2Synonym5.pred(), "hasRelatedSynonym");
        assertEquals(doid2Synonym5.val(), "bedsore");
        assertEquals(doid2Synonym5.xrefs(), Lists.newArrayList());

        DoidBasicPropertyValue doid2BasicPropertyValue1 = doidEntry2.doidMetadata().basicPropertyValues().get(0);
        assertEquals(doid2BasicPropertyValue1.pred(), "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId");
        assertEquals(doid2BasicPropertyValue1.val(), "DOID:8808");

        DoidBasicPropertyValue doid2BasicPropertyValue2 = doidEntry2.doidMetadata().basicPropertyValues().get(1);
        assertEquals(doid2BasicPropertyValue2.pred(), "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId");
        assertEquals(doid2BasicPropertyValue2.val(), "DOID:9129");

        DoidBasicPropertyValue doid2BasicPropertyValue3 = doidEntry2.doidMetadata().basicPropertyValues().get(2);
        assertEquals(doid2BasicPropertyValue3.pred(), "http://www.geneontology.org/formats/oboInOwl#hasOBONamespace");
        assertEquals(doid2BasicPropertyValue3.val(), "disease_ontology");

        DoidBasicPropertyValue doid2BasicPropertyValue4 = doidEntry2.doidMetadata().basicPropertyValues().get(3);
        assertEquals(doid2BasicPropertyValue4.pred(), "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId");
        assertEquals(doid2BasicPropertyValue4.val(), "DOID:9029");

        DoidBasicPropertyValue doid2BasicPropertyValue5 = doidEntry2.doidMetadata().basicPropertyValues().get(4);
        assertEquals(doid2BasicPropertyValue5.pred(), "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId");
        assertEquals(doid2BasicPropertyValue5.val(), "DOID:9002");

    }
}