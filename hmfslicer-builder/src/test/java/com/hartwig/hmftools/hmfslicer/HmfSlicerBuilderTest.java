package com.hartwig.hmftools.hmfslicer;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jooq.Record;
import org.jooq.Result;
import org.junit.Test;

public class HmfSlicerBuilderTest {
    @Test
    public void canGenerateGeneRegionsFromDB() throws IOException, HartwigException, SQLException {
        final Result<Record> queryResults = HmfSlicerBuilderRunner.queryEnsembldb();
        final Set<Object> gene_names = queryResults.stream().map(x -> x.get("gene_name")).collect(Collectors.toSet());
        final List<String> genes = HmfSlicerBuilderRunner.readGeneList();
        assertEquals(genes.size(), gene_names.size());
    }
}
