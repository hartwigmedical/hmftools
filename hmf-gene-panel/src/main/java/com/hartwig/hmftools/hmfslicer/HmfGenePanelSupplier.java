package com.hartwig.hmftools.hmfslicer;

import java.io.IOException;
import java.io.InputStream;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfSlicerFileLoader;

public class HmfGenePanelSupplier {

    public static SortedSetMultimap<String, HmfGenomeRegion> asMap() throws IOException, EmptyFileException {
        final InputStream inputStream = HmfGenePanelSupplier.class.getResourceAsStream("/hmf_gene_panel.tsv");
        return HmfSlicerFileLoader.fromInputStream(inputStream);
    }

    public static List<HmfGenomeRegion> asList() throws IOException, EmptyFileException {
        final List<HmfGenomeRegion> result = Lists.newArrayList(asMap().values());
        Collections.sort(result);
        return result;
    }

}
