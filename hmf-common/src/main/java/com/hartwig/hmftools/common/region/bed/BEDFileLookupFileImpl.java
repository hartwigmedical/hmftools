package com.hartwig.hmftools.common.region.bed;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public class BEDFileLookupFileImpl implements BEDFileLookup {

    private final AbstractFeatureReader<BEDFeature, LineIterator> reader;

    public BEDFileLookupFileImpl(final String bedFile) {
        reader = getFeatureReader(bedFile, new BEDCodec(), true);
    }

    @Override
    public double score(@NotNull final GenomePosition position) throws IOException {
        return reader.query(position.chromosome(), (int) position.position(), (int) position.position())
                .stream()
                .findFirst()
                .map(BEDFeature::getScore)
                .orElse(0f);
    }

    public boolean exists(@NotNull final GenomePosition position) throws IOException {
        return reader.query(position.chromosome(), (int) position.position(), (int) position.position()).stream().findAny().isPresent();
    }

    @Override
    public void close() throws IOException {
        reader.close();
    }
}
