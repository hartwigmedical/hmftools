package com.hartwig.hmftools.sullivan;

import org.jetbrains.annotations.NotNull;

class DefaultNameConverter implements FileNameConverter {

    public String apply(@NotNull String originalFileName) {
        // KODU: Default original name we expect: SAMPLE_FLOWCELL_SX_LANE_Rn_001.fastq.gz
        // KODU: Format we expect to convert to: SAMPLE_FLOWCELL_SX_LANE_001_n.fastq.gz
        // KODU: We also support converting starting from the converted format.
        String nameWithoutExtension = originalFileName.substring(0, originalFileName.indexOf("."));

        String splitRegExp = "_";
        String[] parts = nameWithoutExtension.split(splitRegExp);
        if (parts[5].length() == 1) {
            // KODU: We have converted fastq as input
            return nameWithoutExtension + ".fastq";
        } else {
            String readGroup = parts[4].equals("R1") ? "1" : "2";
            return parts[0] + splitRegExp + parts[1] + splitRegExp + parts[2] + splitRegExp + parts[3] +
                    splitRegExp + parts[5] + splitRegExp + readGroup + ".fastq";
        }
    }
}
