package com.hartwig.hmftools.hmfslicer;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.charset.Charset;
import java.sql.SQLException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.jooq.Record;
import org.jooq.Result;
import org.junit.Test;

public class HmfSlicerBuilderTest {
    @Test
    public void canGenerateSlicer() throws IOException, HartwigException, SQLException {
        final Result<Record> queryResults = HmfSlicerBuilderRunner.queryEnsembldb();
        final int geneNumber = Resources.readLines(Resources.getResource("gene_panel"),
                Charset.defaultCharset()).size();
        assertEquals(geneNumber, queryResults.size());
    }
}
