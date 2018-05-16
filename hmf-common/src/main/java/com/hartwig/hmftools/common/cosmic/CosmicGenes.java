package com.hartwig.hmftools.common.cosmic;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.jetbrains.annotations.NotNull;

public final class CosmicGenes {

    private static final int NAME_COLUMN = 0;
    private static final int DESCRIPTION_COLUMN = 1;
    private static final int ENTREZ_ID_COLUMN = 2;
    private static final int GENOME_LOCATION_COLUMN = 3;
    private static final int CHROMOSOME_BAND_COLUMN = 4;
    private static final int SOMATIC_COLUMN = 5;
    private static final int GERMLINE_COLUMN = 6;
    private static final int SOMATIC_TUMOR_TYPES_COLUMN = 7;
    private static final int GERMLINE_TUMOR_TYPES_COLUMN = 8;
    private static final int CANCER_SYNDROME_COLUMN = 9;
    private static final int TISSUE_TYPE_COLUMN = 10;
    private static final int MOLECULAR_GENETICS_COLUMN = 11;
    private static final int ROLE_COLUMN = 12;
    private static final int MUTATION_TYPES_COLUMN = 13;
    private static final int TRANSLOCATION_PARTNER_COLUMN = 14;
    private static final int OTHER_GERMLINE_MUTATIONS_COLUMN = 15;
    private static final int OTHER_SYNDROME_COLUMN = 16;
    private static final int SYNONYMS_COLUMN = 17;

    private static final String FIELD_SEPARATOR = ",";
    private static final String FIELD_CONCATENATOR = ";";
    private static final String QUOTE = "\"";

    private CosmicGenes() {
    }

    @NotNull
    public static CosmicGeneModel readFromCSV(@NotNull final String pathToCsv) throws IOException {
        final Map<String, CosmicGeneData> dataPerGene = Maps.newHashMap();
        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines) {
            // KODU: Some fields are quoted as they hold the field separator char.
            // Inside these quotes, we replace the field separator by a separate char.
            while (line.contains(QUOTE)) {
                final int startIndex = line.indexOf(QUOTE);
                final int endIndex = line.indexOf(QUOTE, startIndex + 1);
                assert endIndex > startIndex;
                final String field = line.substring(startIndex + 1, endIndex).replace(FIELD_SEPARATOR, FIELD_CONCATENATOR);
                line = line.substring(0, startIndex) + field + line.substring(endIndex + 1);
            }

            final String[] parts = line.split(FIELD_SEPARATOR, 18);
            if (parts.length > 0) {
                final String gene = parts[NAME_COLUMN].trim();
                final CosmicGeneData data = ImmutableCosmicGeneData.of(parts[DESCRIPTION_COLUMN].trim(),
                        parts[ENTREZ_ID_COLUMN].trim(),
                        parts[GENOME_LOCATION_COLUMN].trim(),
                        parts[CHROMOSOME_BAND_COLUMN].trim(),
                        parts[SOMATIC_COLUMN].trim(),
                        parts[GERMLINE_COLUMN].trim(),
                        parts[SOMATIC_TUMOR_TYPES_COLUMN].trim(),
                        parts[GERMLINE_TUMOR_TYPES_COLUMN].trim(),
                        parts[CANCER_SYNDROME_COLUMN].trim(),
                        parts[TISSUE_TYPE_COLUMN].trim(),
                        parts[MOLECULAR_GENETICS_COLUMN].trim(),
                        parts[ROLE_COLUMN].trim(),
                        parts[MUTATION_TYPES_COLUMN].trim(),
                        parts[TRANSLOCATION_PARTNER_COLUMN].trim(),
                        parts[OTHER_GERMLINE_MUTATIONS_COLUMN].trim(),
                        parts[OTHER_SYNDROME_COLUMN].trim());

                dataPerGene.put(gene, data);

                for (final String synonym : parts[SYNONYMS_COLUMN].trim().split(FIELD_CONCATENATOR)) {
                    dataPerGene.put(synonym, data);
                }
            }
        }
        return ImmutableCosmicGeneModel.of(dataPerGene);
    }
}
