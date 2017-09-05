package com.hartwig.hmftools.hmfslicer;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

import org.junit.Test;

public class HmfGenePanelSupplierTest {

    @Test
    public void canLoadGeneRegionsFromFile() throws IOException, EmptyFileException {
        final List<String> genes = HmfGenePanelBuilder.readGeneList();
        final List<HmfGenomeRegion> geneRegions = HmfGenePanelSupplier.asList();
        assertEquals(genes.size(), geneRegions.size());
    }

}
