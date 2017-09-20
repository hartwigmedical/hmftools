import static org.junit.Assert.assertNotNull;

import static junit.framework.TestCase.assertNull;

import java.io.IOException;
import java.nio.file.Paths;

import javax.xml.bind.JAXBException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.bachelor.BachelorSchema;

import nl.hartwigmedicalfoundation.bachelor.Program;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;
import org.xml.sax.SAXException;

public class BachelorTest {

    private final static String TEST_XML = Resources.getResource("valid.xml").getPath();
    private final static String TEST_INVALID_XML = Resources.getResource("missing_namespace.xml").getPath();
    private final static String TEST_DBSNP_XML = Resources.getResource("missing_namespace.xml").getPath();

    @Nullable
    private static Program testFile(final String path) throws JAXBException, SAXException {
        final BachelorSchema schema = BachelorSchema.make();
        return schema.processXML(Paths.get(path));
    }

    @Test
    public void TestValid() throws JAXBException, IOException, SAXException {
        assertNotNull(testFile(TEST_XML));
    }

    @Test
    public void TestMissingNamespace() throws JAXBException, IOException, SAXException {
        assertNull(testFile(TEST_INVALID_XML));
    }

    @Test
    public void TestInvalidDBSNP() throws JAXBException, IOException, SAXException {
        assertNull(testFile(TEST_DBSNP_XML));
    }
}
