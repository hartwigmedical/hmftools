package com.hartwig.hmftools.common.fusions;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Maps;
import com.google.common.collect.Streams;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.math3.util.Pair;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownFusionsModel {
    private static final String FIVE_GENE_COLUMN = "H_gene";
    private static final String THREE_GENE_COLUMN = "T_gene";
    private static final String PROMISCUOUS_GENE_COLUMN = "gene";
    public static final String CGI = "CGI";
    public static final String CIVIC = "CIViC";
    public static final String COSMIC = "COSMIC";
    public static final String ONCOKB = "OncoKB";
    private static final CSVFormat CSV_FORMAT = CSVFormat.DEFAULT.withNullString("NA").withFirstRecordAsHeader();

    //MIVO: fusion pair -> set of sources
    @NotNull
    public abstract Map<Pair<String, String>, Set<String>> fusions();

    //MIVO: promiscuous gene -> set of sources
    @NotNull
    public abstract Map<String, Set<String>> promiscuousFive();

    @NotNull
    public abstract Map<String, Set<String>> promiscuousThree();

    @NotNull
    public static KnownFusionsModel fromInputStreams(@NotNull final InputStream fusionPairsStream,
            @NotNull final InputStream promiscuousFiveStream, @NotNull final InputStream promiscuousThreeStream) throws IOException {
        return ImmutableKnownFusionsModel.of(readFusions(fusionPairsStream),
                readPromiscuous(promiscuousFiveStream),
                readPromiscuous(promiscuousThreeStream));
    }

    @NotNull
    private static Map<Pair<String, String>, Set<String>> readFusions(@NotNull final InputStream stream) throws IOException {
        final Map<Pair<String, String>, Set<String>> fusions = Maps.newHashMap();
        final CSVParser parser = CSVParser.parse(stream, Charset.defaultCharset(), CSV_FORMAT);
        Streams.stream(parser).forEach(record -> {
            final Pair<String, String> fusion = Pair.create(record.get(FIVE_GENE_COLUMN), record.get(THREE_GENE_COLUMN));
            final Set<String> sources = readSources(record);
            fusions.put(fusion, sources);
        });
        return fusions;
    }

    @NotNull
    private static Map<String, Set<String>> readPromiscuous(@NotNull final InputStream stream) throws IOException {
        final CSVParser parser = CSVParser.parse(stream, Charset.defaultCharset(), CSV_FORMAT);
        return Streams.stream(parser)
                .collect(Collectors.toMap(record -> record.get(PROMISCUOUS_GENE_COLUMN), KnownFusionsModel::readSources));
    }

    private boolean exactMatch(@NotNull final String fiveGene, @NotNull final String threeGene) {
        return fusions().containsKey(Pair.create(fiveGene, threeGene));
    }

    public boolean exactMatch(@NotNull final Collection<String> fiveGeneNames, @NotNull final Collection<String> threeGeneNames) {
        return fiveGeneNames.stream().anyMatch(fiveGene -> threeGeneNames.stream().anyMatch(threeGene -> exactMatch(fiveGene, threeGene)));
    }

    public boolean intergenicPromiscuousMatch(@NotNull final Collection<String> fiveGeneNames,
            @NotNull final Collection<String> threeGeneNames) {
        return fiveGeneNames.stream().noneMatch(threeGeneNames::contains) && (
                fiveGeneNames.stream().anyMatch(promiscuousFive()::containsKey) || threeGeneNames.stream()
                        .anyMatch(promiscuousThree()::containsKey));
    }

    public boolean intragenicPromiscuousMatch(@NotNull final Collection<String> fiveGeneNames,
            @NotNull final Collection<String> threeGeneNames) {
        return fiveGeneNames.stream().anyMatch(threeGeneNames::contains) && threeGeneNames.stream()
                .anyMatch(promiscuousThree()::containsKey);
    }

    @NotNull
    public String primarySource(@NotNull final Collection<String> fiveGeneNames, @NotNull final Collection<String> threeGeneNames) {
        final Set<String> sources = fiveGeneNames.stream()
                .flatMap(fiveGene -> threeGeneNames.stream().map(threeGene -> Pair.create(fiveGene, threeGene)))
                .flatMap(fusionPair -> fusions().getOrDefault(fusionPair, Collections.emptySet()).stream())
                .collect(Collectors.toSet());
        return primarySource(sources);
    }

    @NotNull
    private static Set<String> readSources(@NotNull final CSVRecord record) {
        return Stream.of(readSource(record, CGI), readSource(record, CIVIC), readSource(record, COSMIC), readSource(record, ONCOKB))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toSet());
    }

    @NotNull
    private static String primarySource(@NotNull final Set<String> sources) {
        return sources.contains(ONCOKB)
                ? ONCOKB
                : sources.contains(COSMIC) ? COSMIC : sources.contains(CGI) ? CGI : sources.contains(CIVIC) ? CIVIC : "";
    }

    @NotNull
    private static Optional<String> readSource(@NotNull final CSVRecord record, @NotNull final String sourceName) {
        return Optional.ofNullable(record.get(sourceName)).map(value -> sourceName);
    }
}
