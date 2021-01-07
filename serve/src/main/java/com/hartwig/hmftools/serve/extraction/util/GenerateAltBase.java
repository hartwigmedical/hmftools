package com.hartwig.hmftools.serve.extraction.util;

import java.io.File;
import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class GenerateAltBase {

    private GenerateAltBase() {
    }

    private static final Logger LOGGER = LogManager.getLogger(GenerateAltBase.class);

    @NotNull
    public static String createAltOfRefBase(@NotNull String chromosome, Long genomicPosition) throws IOException {
        String extractRefBaseOfPosition = extractRefBaseOfGenomicPosition(chromosome, genomicPosition);

        if (extractRefBaseOfPosition.equals("A")) {
            return  "T";
        } else if (extractRefBaseOfPosition.equals("C")) {
            return  "A";
        } else if (extractRefBaseOfPosition.equals("T")) {
            return  "G";
        } else if (extractRefBaseOfPosition.equals("G")) {
            return  "A";
        } else {
            LOGGER.warn("None alt base can be generate of ref {}", extractRefBaseOfPosition);
            return Strings.EMPTY;
        }
    }

    public static String extractRefBaseOfGenomicPosition(@Nullable String chromosome, Long genomicPosition) throws IOException {
        IndexedFastaSequenceFile fastaSequenceFile =
                new IndexedFastaSequenceFile(new File(System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta"));
        return fastaSequenceFile.getSubsequenceAt(chromosome, genomicPosition, genomicPosition).getBaseString();
    }
}
