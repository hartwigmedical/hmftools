package com.hartwig.hmftools.bachelor;

import static org.junit.Assert.assertNotNull;

import static junit.framework.TestCase.assertNull;

import java.nio.file.Paths;

import com.google.common.io.Resources;

import nl.hartwigmedicalfoundation.bachelor.Program;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;
import org.xml.sax.SAXException;

public class BachelorTest {

    private final static String TEST_XML = Resources.getResource("valid.xml").getPath();
    private final static String TEST_INVALID_XML = Resources.getResource("missing_namespace.xml").getPath();
    private final static String TEST_DBSNP_XML = Resources.getResource("invalid_dbsnp.xml").getPath();

    @Nullable
    private static Program testFile(final String path) throws SAXException {
        final BachelorSchema schema = BachelorSchema.make();
        return schema.processXML(Paths.get(path));
    }

    @Test
    public void TestValid() throws SAXException {
        assertNotNull(testFile(TEST_XML));
    }

    @Test
    public void TestMissingNamespace() throws SAXException {
        assertNull(testFile(TEST_INVALID_XML));
    }

    @Test
    public void TestInvalidDBSNP() throws SAXException {
        assertNull(testFile(TEST_DBSNP_XML));
    }
}
