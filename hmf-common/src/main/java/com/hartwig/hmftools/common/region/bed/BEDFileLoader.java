package com.hartwig.hmftools.common.region.bed;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public enum BEDFileLoader {

    ;
    private static final Logger LOGGER = LogManager.getLogger(BEDFileLoader.class);

    @NotNull
    public static SortedSetMultimap<String, GenomeRegion> fromBedFile(@NotNull String bedFile) throws IOException {
        final SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();

        String prevChromosome = null;
        GenomeRegion prevRegion = null;
        long numberOfBases = 0;
        try (final AbstractFeatureReader<BEDFeature, LineIterator> reader = getFeatureReader(bedFile, new BEDCodec(), false)) {
            for (final BEDFeature bedFeature : reader.iterator()) {
                final String chromosome = bedFeature.getContig();
                final long start = bedFeature.getStart();
                final long end = bedFeature.getEnd();

                if (end < start) {
                    LOGGER.warn("Invalid genome region found in chromosome " + chromosome + ": start=" + start + ", end=" + end);
                } else {
                    final GenomeRegion region = GenomeRegionFactory.create(chromosome, start, end);
                    if (prevRegion != null && chromosome.equals(prevChromosome) && prevRegion.end() >= start) {
                        LOGGER.warn("BED file is not sorted, please fix! Current=" + region + ", Previous=" + prevRegion);
                    } else {
                        regionMap.put(chromosome, region);
                        prevChromosome = chromosome;
                        prevRegion = region;
                        numberOfBases += region.bases();
                    }
                }
            }
        }

        LOGGER.debug("Created slicer from " + bedFile + ": " + regionMap.size() + " regions covering " + numberOfBases + " bases");
        return regionMap;
    }
}
