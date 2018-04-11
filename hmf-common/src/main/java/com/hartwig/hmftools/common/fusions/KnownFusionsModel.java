package com.hartwig.hmftools.common.fusions;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.SetMultimap;
import com.google.common.collect.Streams;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownFusionsModel {
    private static final int FIVE_GENE_COLUMN = 0;
    private static final int THREE_GENE_COLUMN = 1;

    @NotNull
    public abstract SetMultimap<String, String> fusions();

    @NotNull
    public abstract Set<String> promiscuousFive();

    @NotNull
    public abstract Set<String> promiscuousThree();

    public boolean match(@NotNull final String fiveGene, @NotNull final String threeGene) {
        return promiscuousFive().contains(fiveGene) || promiscuousThree().contains(threeGene) || fusions().get(fiveGene)
                .contains(threeGene);
    }

    @NotNull
    public static KnownFusionsModel fromInputStreams(@NotNull final InputStream fusionPairsStream,
            @NotNull final InputStream promiscuousFiveStream, @NotNull final InputStream promiscuousThreeStream) throws IOException {
        return ImmutableKnownFusionsModel.of(readFusions(fusionPairsStream),
                readPromiscuous(promiscuousFiveStream),
                readPromiscuous(promiscuousThreeStream));
    }

    private static SetMultimap<String, String> readFusions(@NotNull final InputStream stream) throws IOException {
        final SetMultimap<String, String> fusionPairs = HashMultimap.create();
        final CSVParser parser = CSVParser.parse(stream, Charset.defaultCharset(), CSVFormat.DEFAULT.withSkipHeaderRecord());
        Streams.stream(parser).forEach(record -> fusionPairs.put(record.get(FIVE_GENE_COLUMN), record.get(THREE_GENE_COLUMN)));
        return fusionPairs;
    }

    private static Set<String> readPromiscuous(@NotNull final InputStream stream) throws IOException {
        final CSVParser parser = CSVParser.parse(stream, Charset.defaultCharset(), CSVFormat.DEFAULT.withSkipHeaderRecord());
        return Streams.stream(parser).map(record -> record.get(0)).collect(Collectors.toSet());
    }
}
