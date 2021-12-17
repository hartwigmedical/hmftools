package com.hartwig.hmftools.common.genome.region;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public final class BEDFileLoader {

    private static final Logger LOGGER = LogManager.getLogger(BEDFileLoader.class);

    private BEDFileLoader() {
    }

    @NotNull
    public static SortedSetMultimap<String, GenomeRegion> fromBedFile(@NotNull String bedFile) throws IOException {
        final SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();

        String prevChromosome = null;
        GenomeRegion prevRegion = null;
        try (final AbstractFeatureReader<BEDFeature, LineIterator> reader = getFeatureReader(bedFile, new BEDCodec(), false)) {
            for (final BEDFeature bedFeature : reader.iterator()) {
                final String chromosome = bedFeature.getContig();
                final int start = bedFeature.getStart();
                final int end = bedFeature.getEnd();

                if (end < start) {
                    LOGGER.warn("Invalid genome region found in chromosome {}: start={}, end={}", chromosome, start, end);
                } else {
                    final GenomeRegion region = GenomeRegions.create(chromosome, start, end);
                    if (prevRegion != null && chromosome.equals(prevChromosome) && prevRegion.end() >= start) {
                        LOGGER.warn("BED file is not sorted, please fix! Current={}, Previous={}", region, prevRegion);
                    } else {
                        regionMap.put(chromosome, region);
                        prevChromosome = chromosome;
                        prevRegion = region;
                    }
                }
            }
        }

        return regionMap;
    }
}
