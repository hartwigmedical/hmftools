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

import org.jetbrains.annotations.NotNull;

public enum HmfGenePanelSupplier {
    ;

    @NotNull
    public static SortedSetMultimap<String, HmfGenomeRegion> defaultMap() throws IOException, EmptyFileException {
        final InputStream inputStream = HmfGenePanelSupplier.class.getResourceAsStream("/hmf_gene_panel.tsv");
        return HmfSlicerFileLoader.fromInputStream(inputStream);
    }

    @NotNull
    public static List<HmfGenomeRegion> defaultList() throws IOException, EmptyFileException {
        return toList(defaultMap());
    }

    @NotNull
    public static List<HmfGenomeRegion> fromFile(@NotNull final String filename) throws IOException, EmptyFileException {
        return toList(HmfSlicerFileLoader.fromFile(filename));
    }

    @NotNull
    private static List<HmfGenomeRegion> toList(@NotNull final SortedSetMultimap<String, HmfGenomeRegion>  map) throws IOException, EmptyFileException {
        final List<HmfGenomeRegion> result = Lists.newArrayList(map.values());
        Collections.sort(result);
        return result;
    }
}
