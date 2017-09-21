package com.hartwig.hmftools.metadata;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Ignore;
import org.junit.Test;

public class MetadataCreatorTest {

    @Test
    @Ignore
    public void works() throws IOException {
        String runsDir = Resources.getResource("runs").getPath();

        MetadataCreator creator = new MetadataCreator(runsDir);
        creator.run();
    }
}