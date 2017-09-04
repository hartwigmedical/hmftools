package com.hartwig.hmftools.hmfslicer;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfSlicerFileLoader;

public class HmfGeneRegionSupplier {

    public static List<HmfGenomeRegion> get() throws IOException, EmptyFileException {
        final InputStream inputStream = HmfGeneRegionSupplier.class.getResourceAsStream("/hmf_gene_panel.tsv");
        return HmfSlicerFileLoader.fromInputStream(inputStream);
    }
}
