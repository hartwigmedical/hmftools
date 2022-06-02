package com.hartwig.hmftools.common.lilac;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.FileReaderUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LilacDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(LilacDataLoader.class);

    private static final String DELIMITER = ",";

    private LilacDataLoader() {
    }

    @NotNull
    public static LilacData load(@NotNull String lilacQcCsv, @NotNull String lilacResultCsv) throws IOException {
        LOGGER.info("Loading LILAC data from {}", new File(lilacQcCsv).getParent());

        String qc = readQC(lilacQcCsv);
        LOGGER.info(" Read QC status '{}' from {}", qc, lilacQcCsv);

        List<LilacRecord> records = readRecords(lilacResultCsv);
        LOGGER.info(" Read {} LILAC records from {}", records.size(), lilacResultCsv);

        return ImmutableLilacData.builder().qc(qc).records(records).build();
    }

    @NotNull
    private static String readQC(@NotNull String lilacQcCsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(lilacQcCsv).toPath());

        Map<String, Integer> fields = FileReaderUtils.createFieldsIndexMap(lines.get(0), DELIMITER);

        String[] values = lines.get(1).split(DELIMITER);

        return values[fields.get("Status")];
    }

    @NotNull
    private static List<LilacRecord> readRecords(@NotNull String lilacResultCsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(lilacResultCsv).toPath());

        Map<String, Integer> fields = FileReaderUtils.createFieldsIndexMap(lines.get(0), DELIMITER);

        List<LilacRecord> records = Lists.newArrayList();
        for (String line : lines.subList(1, lines.size())) {
            String[] values = line.split(DELIMITER);

            int rnaTotal = Integer.parseInt(values[fields.get("RnaTotal")]);

            records.add(ImmutableLilacRecord.builder()
                    .allele(values[fields.get("Allele")])
                    .refFragments(Integer.parseInt(values[fields.get("RefTotal")]))
                    .tumorFragments(Integer.parseInt(values[fields.get("TumorTotal")]))
                    .rnaFragments(rnaTotal > 0 ? rnaTotal : null)
                    .tumorCopyNumber(Double.parseDouble(values[fields.get("TumorCopyNumber")]))
                    .somaticMissense(Double.parseDouble(values[fields.get("SomaticMissense")]))
                    .somaticNonsenseOrFrameshift(Double.parseDouble(values[fields.get("SomaticNonsenseOrFrameshift")]))
                    .somaticSplice(Double.parseDouble(values[fields.get("SomaticSplice")]))
                    .somaticSynonymous(Double.parseDouble(values[fields.get("SomaticSynonymous")]))
                    .somaticInframeIndel(Double.parseDouble(values[fields.get("SomaticInframeIndel")]))
                    .build());
        }

        return records;
    }
}
