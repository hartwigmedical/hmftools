package com.hartwig.hmftools.patientdb.bqr;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class BQRFile {

    private static final String DELIMITER = "\t";

    private BQRFile() {
    }

    @NotNull
    public static List<BQREntry> read(@NotNull String bqrTsv) throws IOException {
        List<BQREntry> entries = Lists.newArrayList();

        List<String> lines = Files.readAllLines(new File(bqrTsv).toPath());

        Map<String, Integer> fields = createFieldsIndexMap(lines.get(0), DELIMITER);

        for (String line : lines.subList(1, lines.size())) {
            String[] values = line.split(DELIMITER, -1);

            entries.add(ImmutableBQREntry.builder()
                    .ref(values[fields.get("ref")])
                    .alt(values[fields.get("alt")])
                    .trinucleotideContext(values[fields.get("trinucleotideContext")])
                    .count(Integer.parseInt(values[fields.get("count")]))
                    .origQuality(Double.parseDouble(values[fields.get("originalQual")]))
                    .recalibratedQuality(Double.parseDouble(values[fields.get("recalibratedQual")]))
                    .build());
        }

        return entries;
    }
}
