package com.hartwig.hmftools.sullivan;

class AncientNameConverter implements FileNameConverter {

    // KODU: Initially (until ~oct2015) we did not use flowcell names in the files.
    public String apply(String originalFileName) {
        // KODU: Default original name we expect: SAMPLE_SX_LANE_Rn_001.fastq.gz
        // KODU: Format we expect to convert to: SAMPLE_SX_LANE_001_n.fastq.gz
        String nameWithoutExtension = originalFileName.substring(0, originalFileName.indexOf("."));

        String splitRegExp = "_";
        String[] parts = nameWithoutExtension.split(splitRegExp);
        String readGroup = parts[3].equals("R1") ? "1" : "2";
        return parts[0] + splitRegExp + parts[1] + splitRegExp + parts[2] +
                splitRegExp + parts[4] + splitRegExp + readGroup + ".fastq";
    }
}
