package com.hartwig.hmftools.serve.sources.vicc;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.serve.actionability.EvidenceLevel;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ActionableEvidenceFactoryTest {

    @Test
    public void canReformatDrugs() {
        assertEquals("Imatinib,Imatinib", ActionableEvidenceFactory.reformatDrugLabels("IMATINIB,IMATINIB"));

        assertNull(ActionableEvidenceFactory.reformatDrugLabels(null));
    }

    @Test
    public void canReformatField() {
        assertEquals("Field", ActionableEvidenceFactory.reformatField("Field"));
        assertEquals("Field", ActionableEvidenceFactory.reformatField("field"));
        assertEquals("Field", ActionableEvidenceFactory.reformatField("FIELD"));

        assertEquals("F", ActionableEvidenceFactory.reformatField("F"));
        assertEquals("F", ActionableEvidenceFactory.reformatField("f"));
        assertEquals("", ActionableEvidenceFactory.reformatField(""));
        assertNull(ActionableEvidenceFactory.reformatField(null));
    }

    @Test
    public void canResolveDirection() {
        assertEquals(EvidenceDirection.RESPONSIVE, ActionableEvidenceFactory.resolveDirection("Responsive"));
        assertEquals(EvidenceDirection.RESPONSIVE, ActionableEvidenceFactory.resolveDirection("Sensitive"));
        assertEquals(EvidenceDirection.RESISTANT, ActionableEvidenceFactory.resolveDirection("Resistant"));

        assertNull(ActionableEvidenceFactory.resolveDirection(null));
        assertNull(ActionableEvidenceFactory.resolveDirection("Conflicting"));
        assertNull(ActionableEvidenceFactory.resolveDirection("This is no direction"));
    }

    @Test
    public void canResolveLevel() {
        assertEquals(EvidenceLevel.A, ActionableEvidenceFactory.resolveLevel("A"));

        assertNull(ActionableEvidenceFactory.resolveLevel(null));
        assertNull(ActionableEvidenceFactory.resolveLevel("XXX"));
    }

    @Test
    public void canExtractDoid() {
        assertEquals("123", ActionableEvidenceFactory.extractDoid("DOID:123"));
        assertEquals(Strings.EMPTY, ActionableEvidenceFactory.extractDoid("SNOMED:123"));
        assertEquals(Strings.EMPTY, ActionableEvidenceFactory.extractDoid("DOID"));
        assertNull(ActionableEvidenceFactory.extractDoid(null));
    }
}